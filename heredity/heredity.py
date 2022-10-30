import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    print("HHIII!!!!")
    probability = 1
    #print(people)
    for person in people:







        # Commpute probability of one gene, two genes, no genes
        # Calculate probability of OR probability if there are parents
        #if person in one_gene:
        #print(people[person]["trait"])
        #print(people[person]["mother"])
        mother = people[person]["mother"]
        father = people[person]["father"]

        #The order that we consider people does not matter
        if people[person]["mother"] == None and people[person]["father"] == None:
            probOne = PROBS["gene"][1]
            probTwo = PROBS["gene"][2]
            probNo = PROBS["gene"][0]
        # OR
        else:
            if people[mother]["name"] in two_genes:
                probGeneOne = 1 - PROBS["mutation"]
            elif people[mother]["name"] in one_gene:
                print("yes")
                probGeneOne = .5
            else:
                probGeneOne = 0 + PROBS["mutation"]
            probNotGeneOne = 1 - probGeneOne
            probGeneOneNotMutate = probGeneOne - (probGeneOne * PROBS["mutation"])
            probNotGeneOneButMutates = PROBS["mutation"] - (PROBS["mutation"] * probNotGeneOne)

            # PROB GENE ONE =
            # PROB IS GENE 1 AND DOES NOT MUTATES
            # OR PROB IS NOT GENE 1 AND DOES MUTATE

            #probGeneOne = probGeneOneNotMutate + probNotGeneOneButMutates

            if people[father]["name"] in two_genes:
                probGeneTwo = 1 - PROBS["mutation"]
            elif people[father]["name"] in one_gene:
                print("yes")
                probGeneTwo = .5
            else:
                probGeneTwo = 0 + PROBS["mutation"]
            probNotGeneTwo = 1 - probGeneTwo
            probGeneTwoNotMutate = probGeneTwo - (probGeneTwo * PROBS["mutation"])
            probNotGeneTwoButMutates = PROBS["mutation"] - (PROBS["mutation"] * probNotGeneTwo) 

            #probGeneTwo = probGeneTwoNotMutate + probNotGeneTwoButMutates
    
            probOne = probGeneOne * (1 - probGeneTwo) + probGeneTwo * (1 - probGeneOne)
            #probOne = (probGeneOne + probGeneTwo) - (probGeneOne * probGeneTwo)
            probTwo = probGeneOne * probGeneTwo
            probNo = (1 - probGeneOne) * (1 - probGeneTwo)
            #probOne = (probGeneOne + probGeneTwo) - (probGeneOne * probGeneTwo)
            #probTwo = probGeneOne * probGeneTwo
            #probNo = (1 - probGeneOne) * (1 - probGeneTwo)

        '''
        # AND
            if people[person]["mother"] == None and people[person]["father"] == None:
                probTwo = PROBS["gene"][2]
            else:
                if people[person]["mother"]["gene"][2] == 1:
                    probGeneOne = 1
                elif people[person]["mother"]["gene"][1] == 1:
                    probGeneOne = .5
                else:
                    probGeneOne = 0
                probNotGeneOne = 1 - probGeneOne
                probGeneOneNotMutate = probGeneOne - (probGeneOne * PROBS["mutation"])
                probNotGeneOneButMutates = PROBS["mutation"] - (PROBS["mutation"] * probNotGeneOne)

                # PROB GENE ONE =
                # PROB IS GENE 1 AND DOES NOT MUTATES
                # OR PROB IS NOT GENE 1 AND DOES MUTATE

                probGeneOne = probGeneOneNotMutate + probNotGeneOneButMutates

                if people[person]["father"]["gene"][2] == 1:
                    probGeneTwo = 1
                elif people[person]["father"]["gene"][1] == 1:
                    probGeneTwo = .5
                else:
                    probGeneTwo = 0
                probNotGeneTwo = 1 - probGeneTwo
                probGeneTwoNotMutate = probGeneTwo - (probGeneTwo * PROBS["mutation"])
                probNotGeneTwoButMutates = PROBS["mutation"] - (PROBS["mutation"] * probNotGeneTwo) 

                probGeneTwo = probGeneTwoNotMutate + probNotGeneTwoButMutates

                probTwo = probGeneOne * probGeneTwo

        # No genes compute
            if people[person]["mother"] == None and people[person]["father"] == None:
                probTwo = PROBS["gene"][0]
            else:
                if people[person]["mother"]["gene"][2] == 1:
                    probGeneOne = 1
                elif people[person]["mother"]["gene"][1] == 1:
                    probGeneOne = .5
                else:
                    probGeneOne = 0
                probNotGeneOne = 1 - probGeneOne
                probGeneOneNotMutate = probGeneOne - (probGeneOne * PROBS["mutation"])
                probNotGeneOneButMutates = PROBS["mutation"] - (PROBS["mutation"] * probNotGeneOne)

                # PROB GENE ONE =
                # PROB IS GENE 1 AND DOES NOT MUTATES
                # OR PROB IS NOT GENE 1 AND DOES MUTATE

                probGeneOne = probGeneOneNotMutate + probNotGeneOneButMutates

                if people[person]["father"]["gene"][2] == 1:
                    probGeneTwo = 1
                elif people[person]["father"]["gene"][1] == 1:
                    probGeneTwo = .5
                else:
                    probGeneTwo = 0
                probNotGeneTwo = 1 - probGeneTwo
                probGeneTwoNotMutate = probGeneTwo - (probGeneTwo * PROBS["mutation"])
                probNotGeneTwoButMutates = PROBS["mutation"] - (PROBS["mutation"] * probNotGeneTwo) 

                probGeneTwo = probGeneTwoNotMutate + probNotGeneTwoButMutates

                probNo = (1 - probGeneOne) * (1 - probGeneTwo)
            
                '''
        # Commpute if that person has trait
        probNoGeneHasTrait = probNo * PROBS["trait"][0][True]
        probOneGeneHasTrait = probOne * PROBS["trait"][1][True]
        probTwoGeneHasTrait = probTwo * PROBS["trait"][2][True]
        probNoGeneNoTrait = probNo * PROBS["trait"][0][False]
        probOneGeneNoTrait = probOne * PROBS["trait"][1][False]
        probTwoGeneNoTrait = probTwo * PROBS["trait"][2][False]
        # The probability of having the trait =
        # the probability of having both genes and having the trait
        # OR the probability of having 1 gene and having the trait
        # OR the probability of having 0 genes and having the trait
        probHasTrait = probNoGeneHasTrait + probOneGeneHasTrait + probTwoGeneHasTrait
        print(probability)
        if person in one_gene and person in have_trait:
            probability *= probOneGeneHasTrait
        if person in one_gene and person not in have_trait:
            probability *= probOneGeneNoTrait
        if person in two_genes and person in have_trait:
            probability *= probTwoGeneHasTrait
        if person in two_genes and person not in have_trait:
            probability *= probTwoGeneNoTrait
        if person not in one_gene and person not in two_genes and person in have_trait:
            probability *= probNoGeneHasTrait
        if person not in one_gene and person not in two_genes and person not in have_trait:
            probability *= probNoGeneNoTrait

        if people[person]["name"] == "Harry" and person in one_gene and person not in have_trait:
            print(f"Person: {person}")
            print(probOne)
            print(f"Chance 1 copies and no trait: {probOneGeneNoTrait}")

        #print(f"Person: {person}")
        #print(f"Chance 1 copies and no trait: {probOneGeneNoTrait}")
        #if person in have_trait:
         #   probability *= probHasTrait
        #else:
        #    probability *= 1 - probHasTrait
    print(f"prob: {probability:.8f}")
    return probability


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    #print(p)
    #print(f"Probabilities: {probabilities}")
    for person in probabilities:
        #if in one gene two gene trait, add p to prob
        if person in one_gene:
            probabilities[person]["gene"][1] += p
        elif person in two_genes:
            probabilities[person]["gene"][2] += p
        else:
            probabilities[person]["gene"][0] += p
        


        if person in have_trait:
            probabilities[person]["trait"][True] += p
        else:
            probabilities[person]["trait"][False] += p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """

    for person in probabilities:
        probTrue = probabilities[person]["trait"][True]
        probFalse = probabilities[person]["trait"][False]

        if probTrue + probFalse != 0:
            normalizeValue = 1 / (probTrue + probFalse)

        probabilities[person]["trait"][True] *= normalizeValue
        probabilities[person]["trait"][False] *= normalizeValue

        probNo = probabilities[person]["gene"][0]
        probOne = probabilities[person]["gene"][1]
        probTwo = probabilities[person]["gene"][2]

        if probNo + probOne + probTwo != 0:
            normalizeValue = 1 / (probNo + probOne + probTwo)

        probabilities[person]["gene"][0] *= normalizeValue
        probabilities[person]["gene"][1] *= normalizeValue
        probabilities[person]["gene"][2] *= normalizeValue


if __name__ == "__main__":
    main()
