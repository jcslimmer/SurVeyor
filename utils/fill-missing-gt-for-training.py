import argparse

parser = argparse.ArgumentParser(
    description="Replace missing genotype values for paired IDs in a file (in-place)."
)
parser.add_argument("file", help="Path to the file to be modified")
args = parser.parse_args()

# Read the file and store each line as a dictionary with keys 'id' and 'gt'
lines = []
with open(args.file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line.strip():
            continue  # skip empty lines
        parts = line.split()
        if len(parts) == 2:
            lines.append({"id": parts[0], "gt": parts[1]})

# Build a mapping from a base ID (ID with _DUP removed if present) to list of indices in 'lines'
base_map = {}
for idx, entry in enumerate(lines):
    id_val = entry["id"]
    if id_val.endswith("_DUP"):
        base_id = id_val[:-4]  # remove the '_DUP' suffix
    else:
        base_id = id_val
    base_map.setdefault(base_id, []).append(idx)

# Define the set of valid genotype values
valid_genotypes = {"0/1", "1/1", "./1"}

# Process each group that has exactly 2 entries (the pair)
for base_id, indices in base_map.items():
    if len(indices) == 2:
        gt_values = [lines[i]["gt"] for i in indices]
        if "./." in gt_values:
            # Determine which line is missing and which has a valid genotype
            missing_index = indices[gt_values.index("./.")]
            other_index = indices[1 - gt_values.index("./.")]
            if lines[other_index]["gt"] in valid_genotypes:
                # Replace the missing genotype
                lines[missing_index]["gt"] = lines[other_index]["gt"]

# Write out the potentially updated lines back to the same file
with open(args.file, "w") as out:
    for entry in lines:
        out.write(f"{entry['id']} {entry['gt']}\n")