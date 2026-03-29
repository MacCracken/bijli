//! AI integration — daimon/hoosh client for bijli.

use serde::{Deserialize, Serialize};

use crate::error::{BijliError, Result};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DaimonConfig {
    pub endpoint: String,
    pub api_key: Option<String>,
}

impl Default for DaimonConfig {
    fn default() -> Self {
        Self {
            endpoint: "http://localhost:8090".into(),
            api_key: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HooshConfig {
    pub endpoint: String,
}

impl Default for HooshConfig {
    fn default() -> Self {
        Self {
            endpoint: "http://localhost:8088".into(),
        }
    }
}

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

    pub fn new(config: DaimonConfig) -> Result<Self> {
        let client = reqwest::Client::builder()
            .timeout(std::time::Duration::from_secs(30))
            .build()
            .map_err(|e| BijliError::InvalidParameter {
                reason: format!("failed to build HTTP client: {e}"),
            })?;
        Ok(Self { config, client })
    }

    pub async fn register_agent(&self) -> Result<String> {
        tracing::debug!(endpoint = %self.config.endpoint, "registering agent with daimon");
        let body = serde_json::json!({
            "name": "bijli",
            "capabilities": ["electric_fields", "magnetic_fields", "maxwell", "charge_dynamics"],
        });
        let resp = self
            .client
            .post(format!("{}/v1/agents/register", self.config.endpoint))
            .json(&body)
            .send()
            .await
            .map_err(|e| BijliError::InvalidParameter {
                reason: format!("registration request failed: {e}"),
            })?;
        let data: serde_json::Value =
            resp.json()
                .await
                .map_err(|e| BijliError::InvalidParameter {
                    reason: format!("invalid registration response: {e}"),
                })?;
        data["agent_id"]
            .as_str()
            .map(String::from)
            .ok_or_else(|| BijliError::InvalidParameter {
                reason: "registration response missing agent_id field".into(),
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
