//! AI integration — daimon/hoosh client for bijli.
//!
//! # Errors
//!
//! Transport, HTTP-status and response-decoding failures are reported as
//! [`BijliError::InvalidParameter`] with a `daimon:`-prefixed reason (the
//! error enum has no dedicated network variant).

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};

/// Maximum accepted size of a daimon response body (1 MiB).
const MAX_RESPONSE_BYTES: usize = 1 << 20;

/// Connection settings for a daimon agent registry.
///
/// `Debug` redacts `api_key`. `Serialize` still writes it in clear text, as
/// it is a configuration round-trip format — keep serialized configs out of
/// logs and world-readable files.
#[derive(Clone, Serialize, Deserialize)]
pub struct DaimonConfig {
    /// Base URL of the daimon service (e.g. `http://localhost:8090`).
    pub endpoint: String,
    /// Optional API key, sent as `Authorization: Bearer <key>`.
    pub api_key: Option<String>,
}

impl std::fmt::Debug for DaimonConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DaimonConfig")
            .field("endpoint", &self.endpoint)
            .field("api_key", &self.api_key.as_ref().map(|_| "<redacted>"))
            .finish()
    }
}

impl Default for DaimonConfig {
    fn default() -> Self {
        Self {
            endpoint: "http://localhost:8090".into(),
            api_key: None,
        }
    }
}

/// Connection settings for a hoosh inference gateway.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HooshConfig {
    /// Base URL of the hoosh service (e.g. `http://localhost:8088`).
    pub endpoint: String,
}

impl Default for HooshConfig {
    fn default() -> Self {
        Self {
            endpoint: "http://localhost:8088".into(),
        }
    }
}

/// HTTP client for registering bijli with daimon.
#[derive(Debug, Clone)]
pub struct DaimonClient {
    config: DaimonConfig,
    client: reqwest::Client,
}

impl DaimonClient {
    /// Returns a reference to the client's configuration.
    #[must_use]
    pub fn config(&self) -> &DaimonConfig {
        &self.config
    }

    /// Build a client with a 30 s request timeout.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if the HTTP client cannot be built.
    pub fn new(config: DaimonConfig) -> Result<Self> {
        tracing::debug!(endpoint = %config.endpoint, "creating daimon client");
        let client = reqwest::Client::builder()
            .timeout(std::time::Duration::from_secs(30))
            .build()
            .map_err(|e| BijliError::InvalidParameter {
                reason: format!("daimon: failed to build HTTP client: {e}"),
            })?;
        Ok(Self { config, client })
    }

    /// Register bijli as an agent and return the assigned `agent_id`.
    ///
    /// Sends `POST {endpoint}/v1/agents/register` with an
    /// `Authorization: Bearer` header when an API key is configured.
    ///
    /// # Errors
    ///
    /// [`BijliError::InvalidParameter`] if the request fails, the server
    /// answers with a non-success (4xx/5xx) status, the body exceeds 1 MiB or
    /// is not JSON, or the JSON has no string `agent_id`.
    pub async fn register_agent(&self) -> Result<String> {
        let endpoint = self.config.endpoint.trim_end_matches('/');
        tracing::debug!(endpoint = %endpoint, "registering agent with daimon");
        let body = serde_json::json!({
            "name": "bijli",
            "capabilities": ["electric_fields", "magnetic_fields", "maxwell", "charge_dynamics"],
        });
        let mut request = self
            .client
            .post(format!("{endpoint}/v1/agents/register"))
            .json(&body);
        if let Some(key) = &self.config.api_key {
            request = request.bearer_auth(key);
        }
        let resp = request
            .send()
            .await
            .map_err(|e| BijliError::InvalidParameter {
                reason: format!("daimon: registration request failed: {e}"),
            })?
            .error_for_status()
            .map_err(|e| BijliError::InvalidParameter {
                reason: format!("daimon: registration rejected: {e}"),
            })?;
        let bytes = read_limited(resp, MAX_RESPONSE_BYTES).await?;
        let data: serde_json::Value =
            serde_json::from_slice(&bytes).map_err(|e| BijliError::InvalidParameter {
                reason: format!("daimon: invalid registration response: {e}"),
            })?;
        let agent_id = data
            .get("agent_id")
            .and_then(serde_json::Value::as_str)
            .map(String::from)
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: "daimon: registration response missing agent_id field".into(),
            })?;
        tracing::debug!(agent_id = %agent_id, "registered with daimon");
        Ok(agent_id)
    }
}

/// Read a response body, failing once it exceeds `limit` bytes.
async fn read_limited(mut resp: reqwest::Response, limit: usize) -> Result<Vec<u8>> {
    let too_large = || BijliError::InvalidParameter {
        reason: format!("daimon: response body exceeds {limit} bytes"),
    };
    if resp.content_length().is_some_and(|len| len > limit as u64) {
        return Err(too_large());
    }
    let mut buf = Vec::new();
    while let Some(chunk) = resp
        .chunk()
        .await
        .map_err(|e| BijliError::InvalidParameter {
            reason: format!("daimon: failed to read response body: {e}"),
        })?
    {
        if buf.len().saturating_add(chunk.len()) > limit {
            return Err(too_large());
        }
        buf.extend_from_slice(&chunk);
    }
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    use tokio::io::{AsyncReadExt, AsyncWriteExt};

    /// Serve one HTTP/1.1 response on localhost; return the base URL and a
    /// handle yielding the raw request head the client sent.
    async fn one_shot_server(
        status_line: &'static str,
        body: String,
    ) -> (String, tokio::task::JoinHandle<String>) {
        let listener = tokio::net::TcpListener::bind("127.0.0.1:0").await.unwrap();
        let addr = listener.local_addr().unwrap();
        let handle = tokio::spawn(async move {
            let (mut sock, _) = listener.accept().await.unwrap();
            let mut req = Vec::new();
            let mut buf = [0u8; 4096];
            // Read until the end of the headers plus the (small) JSON body.
            loop {
                let n = sock.read(&mut buf).await.unwrap();
                if n == 0 {
                    break;
                }
                req.extend_from_slice(&buf[..n]);
                let text = String::from_utf8_lossy(&req);
                if let Some(head_end) = text.find("\r\n\r\n") {
                    let len = text
                        .lines()
                        .find_map(|l| {
                            l.to_ascii_lowercase()
                                .strip_prefix("content-length:")
                                .map(|v| v.trim().parse::<usize>().unwrap())
                        })
                        .unwrap_or(0);
                    if req.len() >= head_end + 4 + len {
                        break;
                    }
                }
            }
            let resp = format!(
                "{status_line}\r\ncontent-type: application/json\r\ncontent-length: {}\r\nconnection: close\r\n\r\n{body}",
                body.len()
            );
            // The client may hang up early (e.g. oversized body); ignore that.
            let _ = sock.write_all(resp.as_bytes()).await;
            sock.shutdown().await.ok();
            String::from_utf8_lossy(&req).into_owned()
        });
        (format!("http://{addr}"), handle)
    }

    #[test]
    fn test_default_config() {
        let c = DaimonConfig::default();
        assert_eq!(c.endpoint, "http://localhost:8090");
        assert!(c.api_key.is_none());
    }

    #[test]
    fn test_hoosh_default() {
        let c = HooshConfig::default();
        assert_eq!(c.endpoint, "http://localhost:8088");
    }

    #[test]
    fn test_daimon_client_new() {
        let config = DaimonConfig::default();
        let client = DaimonClient::new(config);
        assert!(client.is_ok());
    }

    #[test]
    fn test_daimon_config_with_api_key() {
        let c = DaimonConfig {
            endpoint: "https://custom.host:9090".into(),
            api_key: Some("test-key".into()),
        };
        assert_eq!(c.api_key.as_deref(), Some("test-key"));
    }

    #[test]
    fn test_debug_redacts_api_key() {
        let c = DaimonConfig {
            endpoint: "https://custom.host:9090".into(),
            api_key: Some("super-secret-key".into()),
        };
        let shown = format!("{c:?}");
        assert!(!shown.contains("super-secret-key"), "{shown}");
        assert!(shown.contains("redacted") && shown.contains("custom.host"));
        let client = DaimonClient::new(c).unwrap();
        assert!(!format!("{client:?}").contains("super-secret-key"));
    }

    #[test]
    fn test_daimon_config_serde_roundtrip() {
        let c = DaimonConfig {
            endpoint: "http://example.com".into(),
            api_key: Some("key123".into()),
        };
        let json = serde_json::to_string(&c).unwrap();
        let back: DaimonConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(back.endpoint, c.endpoint);
        assert_eq!(back.api_key, c.api_key);
    }

    #[test]
    fn test_hoosh_config_serde_roundtrip() {
        let c = HooshConfig::default();
        let json = serde_json::to_string(&c).unwrap();
        let back: HooshConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(back.endpoint, c.endpoint);
    }

    #[tokio::test]
    async fn test_register_sends_bearer_and_parses_id() {
        let (url, server) =
            one_shot_server("HTTP/1.1 200 OK", r#"{"agent_id":"abc-123"}"#.into()).await;
        let client = DaimonClient::new(DaimonConfig {
            endpoint: format!("{url}/"),
            api_key: Some("s3cret".into()),
        })
        .unwrap();
        let id = client.register_agent().await.unwrap();
        assert_eq!(id, "abc-123");
        let request = server.await.unwrap();
        assert!(
            request.starts_with("POST /v1/agents/register "),
            "{request}"
        );
        assert!(
            request
                .lines()
                .any(|l| l.eq_ignore_ascii_case("authorization: Bearer s3cret")),
            "{request}"
        );
    }

    #[tokio::test]
    async fn test_register_rejects_error_status() {
        // A 500 with a well-formed body used to be accepted.
        let (url, server) = one_shot_server(
            "HTTP/1.1 500 Internal Server Error",
            r#"{"agent_id":"abc-123"}"#.into(),
        )
        .await;
        let client = DaimonClient::new(DaimonConfig {
            endpoint: url,
            api_key: None,
        })
        .unwrap();
        assert!(client.register_agent().await.is_err());
        let request = server.await.unwrap();
        assert!(
            !request.to_ascii_lowercase().contains("authorization:"),
            "{request}"
        );
    }

    #[tokio::test]
    async fn test_register_rejects_oversized_body() {
        let body = format!(
            r#"{{"agent_id":"x","pad":"{}"}}"#,
            "a".repeat(MAX_RESPONSE_BYTES)
        );
        let (url, server) = one_shot_server("HTTP/1.1 200 OK", body).await;
        let client = DaimonClient::new(DaimonConfig {
            endpoint: url,
            api_key: None,
        })
        .unwrap();
        assert!(client.register_agent().await.is_err());
        let _ = server.await;
    }
}
