// created by Zhili
// License MIT
//
extern crate hdf5;
extern crate ndarray;

#[cfg(feature = "blosc")]
use hdf5::filters::blosc_set_nthreads;
use hdf5::{File, H5Type, Result};
use ndarray::{arr2, Array2, Array1, s};
use std::error::Error;
use hdf5::types::FixedAscii;
use std::env;

extern crate flate2;
use flate2::read;
use flate2::write;
use flate2::Compression;
use std::ffi::OsStr;
use std::fs;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use std::time::Instant;
use std::process::exit;

pub fn writer(filename: &str) -> Box<dyn Write> {
    let path = Path::new(filename);
    let file = match fs::File::create(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.description().to_string()),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufWriter::with_capacity(
            128 * 1024,
            write::GzEncoder::new(file, Compression::default()),
        ))
    } else {
        Box::new(BufWriter::with_capacity(128 * 1024, file))
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // read the file
    let head = "##fileformat=VCFv4.2\n##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n";
    let args: Vec<String> = env::args().collect();
    println!("h5ToVCF, convert a special hdf5 format to VCF format. \n Author: Zhili <zhilizheng@outlook.com>");
    println!("\nRun command: {}", args.join(" ")); 
    if args.len() != 3 {
        println!("Usage: \n  h5ToVCF INPUT.hdf5 OUTPUT.vcf.gz");
        exit(1); 
    }

    let filename = &args[1];
    let h5_file = File::open(filename).expect("can't open hdf5 to read");
    
    let outname = &args[2];
    let mut w_file = writer(outname);
    w_file.write_all(head.as_bytes()).expect("can't write to file");

    let ds_geno = h5_file.dataset("genotypes").expect("Error: can't open the genotype, imputed_par_gts");
    let ds_fam = h5_file.dataset("fams").expect("Error: can't open the fam information");
    let ds_id = h5_file.dataset("ids").expect("Error: can't open the id information");
    let ds_bim = h5_file.dataset("bim").expect("Error: can't open the bim information");

    let n_geno = ds_geno.shape()[0];
    let p_geno = ds_geno.shape()[1];
    let m_geno = ds_geno.shape()[2];

    let n_fam = ds_fam.shape()[0];
    let n_id = ds_id.shape()[0];

    let m_bim = ds_bim.shape()[0];
    let c_bim = ds_bim.shape()[1];

    println!("{} individuals, {} markers in hdf5 file", n_geno, m_geno);

    if p_geno < 2 {
        panic!("genotype dimention is incorrect!");
    }

    if n_id != n_fam {
        panic!("IID and FID is not consistent")
    }

    if n_geno != n_fam {
        panic!("sample size in genotype is not consistent with fam information!");
    }

    if m_geno != m_bim {
        panic!("markers in genotype is not consistent with bim information!");
    }

    if c_bim != 5 {
        panic!("Marker information is not full!"); 
    }
    
    // write fam  
    let dt_fam: Vec<String> = ds_fam.read_1d::<FixedAscii<20>>()?.iter().map(|s| format!("{}", String::from(s.as_str())) ).collect();
    let dt_id: Vec<String> = ds_id.read_1d::<FixedAscii<20>>()?.iter().map(|s| format!("{}", String::from(s.as_str())) ).collect();

    let dt_str = dt_fam.iter().zip(&dt_id).map(|(a, b)| format!("{}_{}", a, b)).collect::<Vec<String>>().join("\t");

    let head2 = format!("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", dt_str);
    w_file.write_all(head2.as_bytes())?;

    let mut now = Instant::now();

    // write bim and geno
    for idx in 0..m_geno {
        let cur_marker = idx + 1;
        if cur_marker % 5000 == 0 {
            let new_now = Instant::now();
            let time = new_now.duration_since(now).as_secs();
            println!("{} markers proceeded, elapsed time {} secs, total marker {}", cur_marker, time, m_geno);
            now = Instant::now();
        }
        w_file.write_all(b"\n").expect("Error to write file");
        let dt_bim: Array1<FixedAscii<20>> = ds_bim.read_slice::<FixedAscii<20>, _, _>(s![idx, ..]).expect("read bim error");
        
        let dt_geno: Array1<f32> = ds_geno.read_slice::<f32, _, _>(s![.., 1, idx]).expect("read geno error");
        let dt_str = dt_geno.iter().map(
           |s|{
               if s.is_finite() {
                   let ret = s / 2.0;
                   ret.to_string()
               }else{
                   String::from(".")
               }
               }).collect::<Vec<String>>().join("\t");

       let to_write = format!("{}\t{}\t{}\t{}\t{}\t.\t.\t.\tDS\t{}", dt_bim[0], dt_bim[2], dt_bim[1],  dt_bim[4], dt_bim[3], dt_str);
       w_file.write_all(to_write.as_bytes()).expect("Error to write file");
    }

    println!("Done");

    Ok(())
}



