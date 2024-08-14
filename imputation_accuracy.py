import pandas as pd
import numpy as np


def reader(size, chrom, datatype):
    return pd.read_csv(
        f"/local/storage/jjh276/gastric_imputation/final_datasets/concordance/python/{size}_chr{chrom}_{datatype}_FINAL.traw", sep="\t"
    ).drop(columns=["SNP", "(C)M"])


def sorter(s):
    if s=="CHR":
        return 0
    elif s=="POS":
        return 1
    elif s=="COUNTED":
        return 2
    elif s=="ALT":
        return 3
    else:
        return 4


def compare(chrom, size1, size2, datatype1, datatype2):
    try:
        print(f"\n\nBeginning comparison between {size1} {datatype1} and {size2} {datatype2}")

        # Read in the files as pandas DataFrames
        input1 = reader(size1, chrom, datatype1)
        input2 = reader(size2, chrom, datatype2)
        print("Finished reading both files")

        # Get the set of dogs that are in both DataFrames
        dog_intersection = list(set.intersection(set(input1.columns), set(input2.columns)))
        dog_intersection.sort(key=sorter)  # Acts in place
        print(f"Found {len(dog_intersection)-4} dogs in common")

        # Select only the dogs that are in both DataFrames (and drop POS duplicate rows)
        input1 = input1[dog_intersection].drop_duplicates(subset="POS", keep=False)
        input2 = input2[dog_intersection]
        input2 = input2[input2["POS"].isin(input1["POS"])].drop_duplicates(subset="POS", keep=False)
        input2.reset_index(drop=True, inplace=True)
        input1 = input1[input1["POS"].isin(input2["POS"])]
        input1.reset_index(drop=True, inplace=True)

        # Check to make sure that the resulting DataFrames have the same CHR and POS values
        assert all(input1["CHR"].values == input2["CHR"].values)
        assert all(input1["POS"].values == input2["POS"].values)
        assert all(input1.columns == input2.columns)
        print(f"Found {len(input1['POS'].values)} positions in common")

        differences = pd.DataFrame(index=input1.index, columns=input1.columns[4:], dtype=float)

        print("Processing differences")
        for index, input1_row in input1.iterrows():
            try:
                input2_row = input2.iloc[index, :]
            except:
                print(f"Failed at {index=}")
                print(f"{input1_row=}")
                raise

            # Get the ordering of this row in input1
            input1_order = input1_row[["COUNTED", "ALT"]].values
            input2_order = input2_row[["COUNTED", "ALT"]].values

            if all(input1_order == input2_order):
                # Set this position in `differences` to the difference
                differences.loc[index] = input1_row[4:] - input2_row[4:]
            elif all(input1_order == np.flip(input2_order)):
                # Set this position in `differences` to the reversed difference
                differences.loc[index] = input1_row[4:] - (2 - input2_row[4:])
            else:
                # Warn that the ordering is neither the same nor reversed
                pass #print(f"Ordering for {index} is neither the same nor reversed in `input1` ({input1_order}) and `input2` ({input2_order})")

        # Note how many were not compared
        uncomparable_rows = differences.isnull().any(axis=1).sum()
        print(f"Could not compare {uncomparable_rows} positions")

        # Write the final `differences` dataframe to a TSV file.
        filename = f"differences_chr{chrom}_{size1}{datatype1}_{size2}{datatype2}.tsv"
        print(f"Printing differences to '{filename}'")
        differences.to_csv(filename, sep="\t")

    except Exception as exception:
        # If anything above went wrong, we'll report it here, and print the error (a.k.a "exception"),
        # but otherwise keep going as if nothing went wrong so that other combinations will be processed.
        print("!"*80)
        print("!"*80)
        print(f"Failed when comparing chr{chrom} {size1}{datatype1} {size2}{datatype2}")
        print(f"{exception=}")
        print("!"*80)
        print("!"*80)
        print()


chromosomes = range(1, 38+1)
sizes = ["173k", "220k", "670k"]

for chrom in chromosomes:
    for isize, size1 in enumerate(sizes):
        for size2 in sizes[isize+1:]:
            compare(chrom, size1, size2, "array", "array")
            compare(chrom, size1, size2, "imputed", "imputed")
    for size1 in sizes:
        for size2 in sizes:
            compare(chrom, size1, size2, "array", "imputed")
