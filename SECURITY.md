# Security Policy

## Supported Versions

| Version | Supported |
|---|---|
| 2.0.x (Cyrius) | ✅ |
| 1.x (the Rust crate) | ❌ — retired: Rust 1.2.0 is the port's specification only ([ADR 0003](docs/adr/0003-retire-the-rust-tree.md)) |

## Reporting a Vulnerability

If you discover a security vulnerability in bijli, please report it
responsibly through **GitHub Security Advisories**:

1. Go to the [Security tab](../../security/advisories) of this repository.
2. Click **"Report a vulnerability"**.
3. Fill in the details and submit.

**Do not open a public issue for security vulnerabilities.**

## Response Timeline

| Action | Target |
|---|---|
| Acknowledgement | Within **48 hours** |
| Initial assessment | Within **5 business days** |
| Fix for critical severity | Within **14 days** |
| Fix for high severity | Within **30 days** |
| Fix for moderate/low severity | Next scheduled release |

## Scope

This policy covers bijli's library, `dist/bijli.cyr` and `dist/bijli-ai.cyr`, and its public API. The surfaces that
read data from elsewhere are the Touchstone reader, the JSON readers, the soorat constructors' component names and
the ai client's HTTP exchange; [`docs/development/threat-model.md`](docs/development/threat-model.md) describes what
each defends against and what remains. A wrong numerical result is a bug, not a vulnerability, unless it reads or
writes memory it shouldn't.

Vulnerabilities in the Cyrius toolchain or its stdlib (sakshi, bayan, sandhi and the rest) should be reported to the
Cyrius project, and flagged here if they affect bijli's users.

## Disclosure

We follow coordinated disclosure. Once a fix is released, we will publish a
security advisory crediting the reporter (unless anonymity is requested).
