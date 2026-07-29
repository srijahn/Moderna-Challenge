import time
import RNA

sequences = [
    "GGUGCCGAACAGUAGCACUC",
    "GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG",
    "AUGCGGAUACGGAUCGUAACGCUAGCUAGCUA"
]

print("\n========== Runtime Analysis ==========\n")

for seq in sequences:

    start = time.time()

    structure, mfe = RNA.fold(seq)

    end = time.time()

    runtime = end - start

    print("Length:", len(seq))
    print("Runtime:", round(runtime, 6), "seconds")
    print("MFE:", mfe)
    print()