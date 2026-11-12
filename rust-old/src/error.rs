//! Error types for bijli.

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum BijliError {
    #[error("invalid field dimensions: expected {expected}D, got {got}D")]
    DimensionMismatch { expected: usize, got: usize },

    #[error("charge magnitude is zero — cannot compute field direction")]
    ZeroCharge,

    #[error("singularity: field evaluation at source point (distance = 0)")]
    Singularity,

    #[error("grid resolution too coarse: {cells} cells for domain size {domain_size}")]
    InsufficientResolution { cells: usize, domain_size: f64 },

    #[error("invalid permittivity: {value} (must be > 0)")]
    InvalidPermittivity { value: f64 },

    #[error("invalid permeability: {value} (must be > 0)")]
    InvalidPermeability { value: f64 },

    #[error("division by zero in EM calculation: {context}")]
    DivisionByZero { context: String },

    #[error("invalid parameter: {reason}")]
    InvalidParameter { reason: String },
}

pub type Result<T> = std::result::Result<T, BijliError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dimension_mismatch() {
        let e = BijliError::DimensionMismatch {
            expected: 3,
            got: 2,
        };
        let msg = e.to_string();
        assert!(msg.contains("3D"));
        assert!(msg.contains("2D"));
    }

    #[test]
    fn test_zero_charge() {
        let e = BijliError::ZeroCharge;
        assert!(e.to_string().contains("zero"));
    }

    #[test]
    fn test_singularity() {
        let e = BijliError::Singularity;
        assert!(e.to_string().contains("distance = 0"));
    }

    #[test]
    fn test_insufficient_resolution() {
        let e = BijliError::InsufficientResolution {
            cells: 4,
            domain_size: 100.0,
        };
        assert!(e.to_string().contains("4 cells"));
    }

    #[test]
    fn test_invalid_permittivity() {
        let e = BijliError::InvalidPermittivity { value: -1.0 };
        assert!(e.to_string().contains("-1"));
    }

    #[test]
    fn test_invalid_permeability() {
        let e = BijliError::InvalidPermeability { value: 0.0 };
        assert!(e.to_string().contains("0"));
    }

    #[test]
    fn test_division_by_zero() {
        let e = BijliError::DivisionByZero {
            context: "field magnitude".into(),
        };
        assert!(e.to_string().contains("field magnitude"));
    }

    #[test]
    fn test_invalid_parameter() {
        let e = BijliError::InvalidParameter {
            reason: "negative frequency".into(),
        };
        assert!(e.to_string().contains("negative frequency"));
    }

    #[test]
    fn test_result_alias() {
        let ok: Result<f64> = Ok(1.0);
        assert!(ok.is_ok());
        let err: Result<f64> = Err(BijliError::ZeroCharge);
        assert!(err.is_err());
    }

    #[test]
    fn test_error_is_debug() {
        let e = BijliError::Singularity;
        let debug = format!("{:?}", e);
        assert!(debug.contains("Singularity"));
    }
}
