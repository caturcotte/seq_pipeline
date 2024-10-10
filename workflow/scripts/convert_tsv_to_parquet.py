import pyarrow as pa
import pyarrow.csv as pc
import pyarrow.parquet as pq


if __name__ == "__main__":
    f = snakemake.input[0]
    o = snakemake.output[0]
    chunksize = 10000
    if snakemake.params.short:
        types = pa.schema([
            ('sample', pa.string()),
            ('chromosome', pa.string()),
            ('position', pa.int64()),
            ('reference', pa.string()),
            ('variant', pa.string()),
            ('quality', pa.float64()),
            ('genotype', pa.string()),
            ('depth', pa.int64()),
            ('allele_depth', pa.string())
        ]) 
    else:
        # types = pa.schema([
        #     ('sample', pa.dictionary(pa.int32(), pa.utf8())),
        #     ('chromosome', pa.dictionary(pa.int32(), pa.utf8())),
        #     ('position', pa.int64()),
        #     ('reference', pa.dictionary(pa.int32(), pa.utf8())),
        #     ('variant', pa.dictionary(pa.int32(), pa.utf8())),
        #     ('quality', pa.string()),
        #     ('genotype', pa.dictionary(pa.int32(), pa.utf8())),
        #     ('depth', pa.int32()),
        #     ('allele_depth', pa.string()),
        #     ('phase_set', pa.string())
        # ]) 
        types = pa.schema([
            ('sample', pa.string()),
            ('chromosome', pa.string()),
            ('position', pa.int64()),
            ('reference', pa.string()),
            ('variant', pa.string()),
            ('quality', pa.string()),
            ('genotype', pa.string()),
            ('depth', pa.string()),
            ('allele_depth', pa.string()),
            ('phase_set', pa.string())
        ]) 

    read_options = pc.ReadOptions(
        block_size=30000,
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
        # convert_options=convert_options,
    )
    writer = pq.ParquetWriter(
        o,
        types
    )
    writer.write_table(stream)
    writer.close()
    print(f"{f} finished.")
