from math import factorial

dominant = int(input())
heterozygous = int(input())
recessive = int(input())
total = dominant + heterozygous + recessive
total_pairs = factorial(total) // (2*factorial(total-2))
dominant_pairings = dominant*(recessive+heterozygous) +factorial(dominant) // (2*factorial(dominant-2))
het_rec_pairs = heterozygous*recessive
het_het_pairs = factorial(heterozygous) // (2*factorial(heterozygous-2))
dominant_p = dominant_pairings/total_pairs
het_rec_p = (het_rec_pairs * 0.5) /total_pairs
het_het_p = (het_het_pairs * 0.75) /total_pairs
print(dominant_p+het_het_p+het_rec_p)
