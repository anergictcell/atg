use std::fs::File;

use s3reader::{S3ObjectUri, S3Reader};

use atglib::utils::errors::AtgError;

// There will be only a single instance of this enum
// so we can allow a large variant
#[allow(clippy::large_enum_variant)]
/// ReadSeekWrapper is an enum to allow dynamic assignment of either File or S3 Readers
/// to be used in the Reader objects of Atglib.
pub enum ReadSeekWrapper {
    File(File, String),
    S3(S3Reader, String),
}

impl ReadSeekWrapper {
    pub fn from_filename(filename: &str) -> Result<Self, AtgError> {
        if filename.starts_with("s3://") {
            let uri = S3ObjectUri::new(filename).map_err(AtgError::new)?;
            let s3obj = S3Reader::open(uri).map_err(AtgError::new)?;
            Ok(Self::S3(s3obj, filename.to_string()))
        } else {
            Ok(Self::File(File::open(filename)?, filename.to_string()))
        }
    }

    pub fn from_cli_arg(filename: &Option<&str>) -> Result<ReadSeekWrapper, AtgError> {
        if let Some(filename) = filename {
            Ok(ReadSeekWrapper::from_filename(filename)?)
        } else {
            Err(AtgError::new("No file specified"))
        }
    }

    pub fn filename(&self) -> &str {
        match self {
            ReadSeekWrapper::File(_, fname) => fname,
            ReadSeekWrapper::S3(_, fname) => fname,
        }
    }
}

// forward all custom implementations straight to the actual reader
impl std::io::Read for ReadSeekWrapper {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize, std::io::Error> {
        match self {
            ReadSeekWrapper::S3(r, _) => r.read(buf),
            ReadSeekWrapper::File(r, _) => r.read(buf),
        }
    }

    fn read_to_end(&mut self, buf: &mut Vec<u8>) -> Result<usize, std::io::Error> {
        match self {
            ReadSeekWrapper::S3(r, _) => r.read_to_end(buf),
            ReadSeekWrapper::File(r, _) => r.read_to_end(buf),
        }
    }

    fn read_to_string(&mut self, buf: &mut String) -> Result<usize, std::io::Error> {
        match self {
            ReadSeekWrapper::S3(r, _) => r.read_to_string(buf),
            ReadSeekWrapper::File(r, _) => r.read_to_string(buf),
        }
    }
}

impl std::io::Seek for ReadSeekWrapper {
    fn seek(&mut self, pos: std::io::SeekFrom) -> Result<u64, std::io::Error> {
        match self {
            ReadSeekWrapper::S3(r, _) => r.seek(pos),
            ReadSeekWrapper::File(r, _) => r.seek(pos),
        }
    }
}
