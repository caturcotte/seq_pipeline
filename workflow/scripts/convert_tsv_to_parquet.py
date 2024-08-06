#!/usr/bin/env python3

import pyarrow as pa
import pyarrow.csv as pc
import pyarrow.parquet as pq

if __name__ == "__main__":
    f = snakemake.input[0]
    o = snakemake.output[0]
    chunksize = 10000
    types = pa.schema([
        ('sample', pa.dictionary(pa.int32(), pa.utf8())),
        ('chromosome', pa.dictionary(pa.int32(), pa.utf8())),
        ('position', pa.int64()),
        ('reference', pa.dictionary(pa.int32(), pa.utf8())),
        ('variant', pa.dictionary(pa.int32(), pa.utf8())),
        ('quality', pa.float32()),
        ('genotype', pa.dictionary(pa.int32(), pa.utf8())),
        ('depth', pa.int32()),
        ('allele_depth', pa.string()),
    ]) 
    read_options = pc.ReadOptions(
        block_size=10000,
    )
    parse_options = pc.ParseOptions(
        delimiter='\t'
    )
    convert_options=pc.ConvertOptions(
        column_types=types,
    )
    # for f in all_files:
    print(f"converting {f} to {o}...")
    stream = pc.read_csv(
        f,
        read_options=read_options,
        parse_options=parse_options,
        convert_options=convert_options,
    )
    writer = pq.ParquetWriter(
        o,
        types
    )
    writer.write_table(stream)
    writer.close()
    print(f"{f} finished.")
