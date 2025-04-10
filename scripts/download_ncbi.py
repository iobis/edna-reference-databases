from Bio import Entrez, SeqIO
import time


def batch(iterable, size):
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]


Entrez.email = "pieterprovoost@gmail.com"
output_file = "sequences.fasta"
query = '(12S[All Fields] OR 16S[All Fields]) AND ribosomal[All Fields] AND ("50"[SLEN] : "50000"[SLEN])'
batch_size = 1000
retmax = 1000000000


open(output_file, "w").close()
search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)
search_results = Entrez.read(search_handle)
search_handle.close()
seq_ids = search_results["IdList"]
print(f"Found {len(seq_ids)} sequences.")

total_batches = sum(1 for _ in batch(seq_ids, batch_size))
start_time = time.time()

for i, batch_items in enumerate(batch(seq_ids, size=batch_size), start=1):
    batch_start = time.time()

    while True:
        try:
            fetch_handle = Entrez.efetch(db="nucleotide", id=",".join(batch_items), rettype="fasta", retmode="text")
            with open(output_file, "a") as out_handle:
                records = SeqIO.parse(fetch_handle, "fasta")
                count = SeqIO.write(records, out_handle, "fasta")
            fetch_handle.close()
            break
        except Exception as e:
            print(e)
            time.sleep(10)

    batch_time = time.time() - batch_start
    elapsed_time = time.time() - start_time
    avg_time_per_batch = elapsed_time / i
    remaining_batches = total_batches - i
    estimated_time_remaining = remaining_batches * avg_time_per_batch

    print(f"Batch {i}/{total_batches}: Saved {count} sequences in {batch_time:.2f} sec.")
    print(f"Estimated time remaining: {estimated_time_remaining:.2f} sec.")
