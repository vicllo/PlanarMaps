import random
from MapPermutation import MapPermutation
from LabelledMap import LabelledMap


def generateRandomBitstring(n):
    L = 4 * n - 2
    weight = n - 1
    bits = [0] * L
    one_positions = random.sample(range(L), weight)
    for pos in one_positions:
        bits[pos] = 1
    return bits


def checkPrefixCondition(bits):
    current_sum = 0
    for j in range(len(bits)-1):
        bit = bits[j]
        current_sum += 3 if bit == 1 else -1
        if current_sum <= -2:
            return False
    return True


def cyclicShift(bits, shift):
    L = len(bits)
    return bits[shift:] + bits[:shift]


"""
def generate_valid_codeword_uniform(n):
    L = 4 * n - 2
    ones_total = n - 1
    zeros_total = 3 * n - 1  # puisque 4n-2 - (n-1) = 3n-1

    # dp(s, r1, r0) = nombre de complétions valides à partir de l'état (s, r1, r0)
    # Condition : tant qu'il reste des coups (r1 + r0 > 0), on exige s >= -1.
    # Au final, quand r1 == r0 == 0, on exige s == -2.
    @lru_cache(maxsize=None)
    def dp(s, r1, r0):
        if r1 == 0 and r0 == 0:
            return 1 if s == -2 else 0
        total = 0
        # Essai de placer un 1 (si disponible) : le pas est +3.
        if r1 > 0:
            new_s = s + 3
            # Si ce n'est pas le dernier coup, il faut que new_s >= -1.
            if (r1 + r0 - 1) > 0:
                if new_s >= -1:
                    total += dp(new_s, r1 - 1, r0)
            else:
                # Si c'est le dernier coup, new_s doit être exactement -2.
                if new_s == -2:
                    total += dp(new_s, r1 - 1, r0)
        # Essai de placer un 0 (si disponible) : le pas est -1.
        if r0 > 0:
            new_s = s - 1
            if (r1 + r0 - 1) > 0:
                if new_s >= -1:
                    total += dp(new_s, r1, r0 - 1)
            else:
                if new_s == -2:
                    total += dp(new_s, r1, r0 - 1)
        return total

    total_count = dp(0, ones_total, zeros_total)
    if total_count == 0:
        raise Exception("Aucune séquence valide n'existe (problème de paramètres)")

    # Maintenant, on génère la séquence de manière récursive.
    seq = []
    s = 0
    r1 = ones_total
    r0 = zeros_total
    for i in range(L):
        count1 = 0
        count0 = 0
        # Option "1" (descente)
        if r1 > 0:
            new_s = s + 3
            if (r1 + r0 - 1) > 0:
                if new_s >= -1:
                    count1 = dp(new_s, r1 - 1, r0)
            else:
                if new_s == -2:
                    count1 = dp(new_s, r1 - 1, r0)
        # Option "0" (remontée)
        if r0 > 0:
            new_s = s - 1
            if (r1 + r0 - 1) > 0:
                if new_s >= -1:
                    count0 = dp(new_s, r1, r0 - 1)
            else:
                if new_s == -2:
                    count0 = dp(new_s, r1, r0 - 1)
        total_choice = count1 + count0
        if total_choice == 0:
            raise Exception("Erreur: aucun mouvement valide possible.")
        r = random.randrange(total_choice)
        if r < count1:
            seq.append(1)
            s += 3
            r1 -= 1
        else:
            seq.append(0)
            s -= 1
            r0 -= 1
    return seq


def test_uniform_sampling(n, num_samples):
    freq = defaultdict(int)
    for _ in range(num_samples):
        cw = tuple(generate_valid_codeword_uniform(n))
        freq[cw] += 1
    return freq

# Test pour n = 3 sur 1 000 000 d'échantillons
n = 3
num_samples = 10**6
frequencies = test_uniform_sampling(n, num_samples)

print("Nombre de codewords distincts pour n = {}: {}".format(n, len(frequencies)))
for codeword, count in sorted(frequencies.items()):
    # Affichage du codeword sous forme de chaîne
    code_str = ''.join(str(bit) for bit in codeword)
    print("Codeword: {} | Fréquence: {} (Relative: {:.5f}%)".format(code_str, count, count/num_samples*100))"""


def generateValidCodeword(n):
    L = 4 * n - 2
    while True:
        b = generateRandomBitstring(n)
        candidate_shifts = [i for i in range(L)]
        valid_shifts = []
        for shift in candidate_shifts:
            if checkPrefixCondition(cyclicShift(b, shift)):
                valid_shifts.append(shift)
        if valid_shifts:
            chosen = random.choice(valid_shifts)
            return cyclicShift(b, chosen)


# Example usage for n = 3:
"""

n = 3
codeword = generateValidCodeword(n)
print("Valid bitstring (length {}):".format(4*n-2), codeword)
L = [1, 0, 0, 1, 0, 0, 0, 0, 0, 0]
count = 0
for i in range(1000000):
    G = generateValidCodeword(n)
    if not checkPrefixCondition(G) :
        print("stop")
        break
    if L == G:
        count +=1
print(count/1000000)"""


def transitiveCouplePermutation(sigma, alpha):
    """
    Check that sigma and alpha act transitively
    """
    assert alpha.size() == sigma.size()
    size = sigma.size()
    seen = [False] * (size + 1)
    seen[0] = seen[1] = True
    # Half-edges are numbered from 1 to size, included

    todo = [1]
    while todo:
        i = todo.pop()
        if not seen[alpha(i)]:
            todo.append(alpha(i))
            seen[alpha(i)] = True
        if not seen[sigma(i)]:
            todo.append(sigma(i))
            seen[sigma(i)] = True
    return False not in seen


def bitToTree(n):
    b = generateValidCodeword(n)
    alphaCycle = [(i, i+1) for i in range(1, 6*n-1, 2)]
    sigmaCycle = []

    p = [1]  # pile contenant les points fixes non couplés
    h = 0  # hauteur dans l'arbre

    H = [0]  # pile des hauteurs associés aux elements de p
    D = [0]*(6*n-1)  # hauteur de toutes les demi-edges
    Q = [[] for _ in range(n)]
    Q[0].append(2)
    D[1] = 0
    D[2] = 0
    Q.append([1])
    j = 4
    i = 0
    while i < len(b) and j <= 6 * n - 2:
        if b[i] == 1:
            D[j] = h+1
            D[j-1] = h
            Q[h].append(j-1)
            Q[h+1].append(j)
            j += 2
            h += 1
            i += 1
            continue
        # b[i] == 0
        if H and H[-1] == h:
            D[j] = h
            D[j-1] = h
            Q[h].append(j-1)
            Q.append([j])
            p.pop()
            H.pop()
            if i+1 < len(b) and b[i+1] == 0:
                if Q[h]:
                    sigmaCycle.append(tuple(Q[h]))
                    Q[h] = []
                h -= 1
            j += 2
            i += 1
            continue
        D[j] = h
        D[j-1] = h
        Q[h].append(j-1)
        p.append(j)
        Q.append([j])
        H.append(h)
        j += 2
        i += 1
    # pour tout h si Q[h] non vide créer un cycle

    for level in range(len(Q)):
        if Q[level]:
            sigmaCycle.append(tuple(Q[level]))
    sigma = MapPermutation(sigmaCycle)
    alpha = MapPermutation(alphaCycle)
    sigma.pretty_print()
    alpha.pretty_print()
    return
    print(sigmaCycle)
    print(alphaCycle)
    return LabelledMap(alpha=alpha, sigma=sigma)


def bitToTreeDebug(b):
    print(len(b))
    n = (len(b)+2)//4
    print(n)

    alphaCycle = [(i, i+1) for i in range(1, 6*n-1, 2)]
    sigmaCycle = []

    p = [1]  # pile contenant les points fixes non couplés
    h = 0  # hauteur dans l'arbre

    H = [0]  # pile des hauteurs associés aux elements de p
    D = [0]*(6*n-1)  # hauteur de toutes les demi-edges
    Q = [[] for _ in range(n)]
    Q[0].append(2)
    D[1] = 0
    D[2] = 0
    Q.append([1])
    j = 4
    i = 0
    while i < len(b) and j <= 6 * n - 2:
        if b[i] == 1:
            D[j] = h+1
            D[j-1] = h
            Q[h].append(j-1)
            Q[h+1].append(j)
            j += 2
            h += 1
            i += 1
            continue
        # b[i] == 0
        if H and H[-1] == h:
            D[j] = h
            D[j-1] = h
            Q[h].append(j-1)
            Q.append([j])
            p.pop()
            H.pop()
            if i+1 < len(b) and b[i+1] == 0:
                if Q[h]:
                    sigmaCycle.append(tuple(Q[h]))
                    Q[h] = []
                h -= 1
            j += 2
            i += 1
            continue
        D[j] = h
        D[j-1] = h
        Q[h].append(j-1)
        p.append(j)
        Q.append([j])
        H.append(h)
        j += 2
        i += 1
    # pour tout h si Q[h] non vide créer un cycle

    for level in range(len(Q)):
        if Q[level]:
            sigmaCycle.append(tuple(Q[level]))
    sigma = MapPermutation(sigmaCycle)
    alpha = MapPermutation(alphaCycle)
    return sigma, alpha
