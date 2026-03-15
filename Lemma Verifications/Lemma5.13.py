import time
import numpy as np
from numba import njit

MAX_N = 8_000_000_000
START_N = 4_810_000_000  # already verified up to and including this
BLOCK = 10_000_000


def _miller_rabin_witness(a, s, d, n):
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return False
    for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1:
            return False
    return True


def is_prime_u64(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p: return True
        if n % p == 0: return False
    d, s = n - 1, 0
    while d % 2 == 0:
        s += 1
        d //= 2
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0: continue
        if _miller_rabin_witness(a, s, d, n): return False
    return True


def squarefree_table(limit):
    sf = np.ones(limit + 1, dtype=np.uint8)
    sf[0] = 0
    p = 2
    while p * p <= limit:
        sf[p*p::p*p] = 0
        p += 1
    return sf


@njit(cache=True)
def gcd_numba(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=True)
def verify_interval_kernel(n_start, n_end, primes_desc, squarefree, all_n, excluded):
    out = np.empty(n_end - n_start + 1, dtype=np.int64)
    bad = 0
    for n in range(n_start, n_end + 1):
        if all_n == 0 and n % 2 == 0: continue
        g = 0
        success = False
        for idx in range(primes_desc.size):
            p = primes_desc[idx]
            d = n - p
            if d <= 0 or d >= squarefree.size: continue
            if squarefree[d] and excluded[d] == 0:
                g = d if g == 0 else gcd_numba(g, d)
                if g <= 2:
                    success = True
                    break
        if not success:
            out[bad] = n
            bad += 1
    return out[:bad]


def largest_primes_below(hi, count=100):
    found = []
    x = hi - 1 if hi % 2 == 0 else hi
    while len(found) < count:
        if is_prime_u64(x):
            found.append(x)
        x -= 2
    return np.array(found, dtype=np.int64)


def main():
    t0 = time.time()
    sf = squarefree_table(2 * BLOCK)
    excluded = np.zeros(2 * BLOCK + 1, dtype=np.uint8)
    exceptions = []

    # START_N = 481 * BLOCK, so the first block has prev_hi = START_N, cur_start = START_N + 1
    a_start = START_N // BLOCK - 1
    max_a = (MAX_N - 1) // BLOCK - 1
    for a in range(a_start, max_a + 1):
        prev_lo = a * BLOCK + 1
        prev_hi = (a + 1) * BLOCK
        cur_start = prev_hi + 1
        cur_end = min((a + 2) * BLOCK, MAX_N)
        if cur_start > cur_end:
            break
        pool = largest_primes_below(prev_hi, 100)
        bad = verify_interval_kernel(cur_start, cur_end, pool, sf, 0, excluded)
        exceptions.extend(int(x) for x in bad)
        print(f"  done block a={a}, up to {cur_end:,}  [{time.time()-t0:.0f}s]")

    exceptions = sorted(set(exceptions))
    print(f"Checked up to {MAX_N:,} in {time.time()-t0:.1f}s")
    if exceptions:
        print(f"Exceptions ({len(exceptions)}): {exceptions}")
    else:
        print("No exceptions found.")


if __name__ == "__main__":
    main()

"""
Output:

done block a=480, up to 4,820,000,000  [6s]
done block a=481, up to 4,830,000,000  [6s]
done block a=482, up to 4,840,000,000  [6s]
done block a=483, up to 4,850,000,000  [6s]
done block a=484, up to 4,860,000,000  [7s]
done block a=485, up to 4,870,000,000  [7s]
done block a=486, up to 4,880,000,000  [7s]
done block a=487, up to 4,890,000,000  [7s]
done block a=488, up to 4,900,000,000  [8s]
done block a=489, up to 4,910,000,000  [8s]
done block a=490, up to 4,920,000,000  [8s]
done block a=491, up to 4,930,000,000  [8s]
done block a=492, up to 4,940,000,000  [9s]
done block a=493, up to 4,950,000,000  [9s]
done block a=494, up to 4,960,000,000  [9s]
done block a=495, up to 4,970,000,000  [9s]
done block a=496, up to 4,980,000,000  [9s]
done block a=497, up to 4,990,000,000  [10s]
done block a=498, up to 5,000,000,000  [10s]
done block a=499, up to 5,010,000,000  [10s]
done block a=500, up to 5,020,000,000  [10s]
done block a=501, up to 5,030,000,000  [11s]
done block a=502, up to 5,040,000,000  [11s]
done block a=503, up to 5,050,000,000  [11s]
done block a=504, up to 5,060,000,000  [11s]
done block a=505, up to 5,070,000,000  [12s]
done block a=506, up to 5,080,000,000  [12s]
done block a=507, up to 5,090,000,000  [12s]
done block a=508, up to 5,100,000,000  [12s]
done block a=509, up to 5,110,000,000  [13s]
done block a=510, up to 5,120,000,000  [13s]
done block a=511, up to 5,130,000,000  [13s]
done block a=512, up to 5,140,000,000  [13s]
done block a=513, up to 5,150,000,000  [13s]
done block a=514, up to 5,160,000,000  [14s]
done block a=515, up to 5,170,000,000  [14s]
done block a=516, up to 5,180,000,000  [14s]
done block a=517, up to 5,190,000,000  [14s]
done block a=518, up to 5,200,000,000  [15s]
done block a=519, up to 5,210,000,000  [15s]
done block a=520, up to 5,220,000,000  [15s]
done block a=521, up to 5,230,000,000  [15s]
done block a=522, up to 5,240,000,000  [16s]
done block a=523, up to 5,250,000,000  [16s]
done block a=524, up to 5,260,000,000  [16s]
done block a=525, up to 5,270,000,000  [16s]
done block a=526, up to 5,280,000,000  [16s]
done block a=527, up to 5,290,000,000  [17s]
done block a=528, up to 5,300,000,000  [17s]
done block a=529, up to 5,310,000,000  [17s]
done block a=530, up to 5,320,000,000  [17s]
done block a=531, up to 5,330,000,000  [18s]
done block a=532, up to 5,340,000,000  [18s]
done block a=533, up to 5,350,000,000  [18s]
done block a=534, up to 5,360,000,000  [18s]
done block a=535, up to 5,370,000,000  [19s]
done block a=536, up to 5,380,000,000  [19s]
done block a=537, up to 5,390,000,000  [19s]
done block a=538, up to 5,400,000,000  [19s]
done block a=539, up to 5,410,000,000  [19s]
done block a=540, up to 5,420,000,000  [20s]
done block a=541, up to 5,430,000,000  [20s]
done block a=542, up to 5,440,000,000  [20s]
done block a=543, up to 5,450,000,000  [20s]
done block a=544, up to 5,460,000,000  [21s]
done block a=545, up to 5,470,000,000  [21s]
done block a=546, up to 5,480,000,000  [21s]
done block a=547, up to 5,490,000,000  [21s]
done block a=548, up to 5,500,000,000  [22s]
done block a=549, up to 5,510,000,000  [22s]
done block a=550, up to 5,520,000,000  [22s]
done block a=551, up to 5,530,000,000  [22s]
done block a=552, up to 5,540,000,000  [22s]
done block a=553, up to 5,550,000,000  [23s]
done block a=554, up to 5,560,000,000  [23s]
done block a=555, up to 5,570,000,000  [23s]
done block a=556, up to 5,580,000,000  [23s]
done block a=557, up to 5,590,000,000  [23s]
done block a=558, up to 5,600,000,000  [24s]
done block a=559, up to 5,610,000,000  [24s]
done block a=560, up to 5,620,000,000  [24s]
done block a=561, up to 5,630,000,000  [24s]
done block a=562, up to 5,640,000,000  [25s]
done block a=563, up to 5,650,000,000  [25s]
done block a=564, up to 5,660,000,000  [25s]
done block a=565, up to 5,670,000,000  [25s]
done block a=566, up to 5,680,000,000  [26s]
done block a=567, up to 5,690,000,000  [26s]
done block a=568, up to 5,700,000,000  [26s]
done block a=569, up to 5,710,000,000  [26s]
done block a=570, up to 5,720,000,000  [26s]
done block a=571, up to 5,730,000,000  [27s]
done block a=572, up to 5,740,000,000  [27s]
done block a=573, up to 5,750,000,000  [27s]
done block a=574, up to 5,760,000,000  [27s]
done block a=575, up to 5,770,000,000  [28s]
done block a=576, up to 5,780,000,000  [28s]
done block a=577, up to 5,790,000,000  [28s]
done block a=578, up to 5,800,000,000  [28s]
done block a=579, up to 5,810,000,000  [29s]
done block a=580, up to 5,820,000,000  [29s]
done block a=581, up to 5,830,000,000  [29s]
done block a=582, up to 5,840,000,000  [29s]
done block a=583, up to 5,850,000,000  [30s]
done block a=584, up to 5,860,000,000  [30s]
done block a=585, up to 5,870,000,000  [30s]
done block a=586, up to 5,880,000,000  [30s]
done block a=587, up to 5,890,000,000  [31s]
done block a=588, up to 5,900,000,000  [31s]
done block a=589, up to 5,910,000,000  [31s]
done block a=590, up to 5,920,000,000  [31s]
done block a=591, up to 5,930,000,000  [32s]
done block a=592, up to 5,940,000,000  [32s]
done block a=593, up to 5,950,000,000  [32s]
done block a=594, up to 5,960,000,000  [32s]
done block a=595, up to 5,970,000,000  [32s]
done block a=596, up to 5,980,000,000  [33s]
done block a=597, up to 5,990,000,000  [33s]
done block a=598, up to 6,000,000,000  [33s]
done block a=599, up to 6,010,000,000  [33s]
done block a=600, up to 6,020,000,000  [34s]
done block a=601, up to 6,030,000,000  [34s]
done block a=602, up to 6,040,000,000  [34s]
done block a=603, up to 6,050,000,000  [34s]
done block a=604, up to 6,060,000,000  [35s]
done block a=605, up to 6,070,000,000  [35s]
done block a=606, up to 6,080,000,000  [35s]
done block a=607, up to 6,090,000,000  [35s]
done block a=608, up to 6,100,000,000  [35s]
done block a=609, up to 6,110,000,000  [36s]
done block a=610, up to 6,120,000,000  [36s]
done block a=611, up to 6,130,000,000  [36s]
done block a=612, up to 6,140,000,000  [36s]
done block a=613, up to 6,150,000,000  [37s]
done block a=614, up to 6,160,000,000  [37s]
done block a=615, up to 6,170,000,000  [37s]
done block a=616, up to 6,180,000,000  [37s]
done block a=617, up to 6,190,000,000  [38s]
done block a=618, up to 6,200,000,000  [38s]
done block a=619, up to 6,210,000,000  [38s]
done block a=620, up to 6,220,000,000  [38s]
done block a=621, up to 6,230,000,000  [39s]
done block a=622, up to 6,240,000,000  [39s]
done block a=623, up to 6,250,000,000  [39s]
done block a=624, up to 6,260,000,000  [39s]
done block a=625, up to 6,270,000,000  [39s]
done block a=626, up to 6,280,000,000  [40s]
done block a=627, up to 6,290,000,000  [40s]
done block a=628, up to 6,300,000,000  [40s]
done block a=629, up to 6,310,000,000  [40s]
done block a=630, up to 6,320,000,000  [41s]
done block a=631, up to 6,330,000,000  [41s]
done block a=632, up to 6,340,000,000  [41s]
done block a=633, up to 6,350,000,000  [41s]
done block a=634, up to 6,360,000,000  [42s]
done block a=635, up to 6,370,000,000  [42s]
done block a=636, up to 6,380,000,000  [42s]
done block a=637, up to 6,390,000,000  [42s]
done block a=638, up to 6,400,000,000  [42s]
done block a=639, up to 6,410,000,000  [43s]
done block a=640, up to 6,420,000,000  [43s]
done block a=641, up to 6,430,000,000  [43s]
done block a=642, up to 6,440,000,000  [43s]
done block a=643, up to 6,450,000,000  [44s]
done block a=644, up to 6,460,000,000  [44s]
done block a=645, up to 6,470,000,000  [44s]
done block a=646, up to 6,480,000,000  [44s]
done block a=647, up to 6,490,000,000  [44s]
done block a=648, up to 6,500,000,000  [45s]
done block a=649, up to 6,510,000,000  [45s]
done block a=650, up to 6,520,000,000  [45s]
done block a=651, up to 6,530,000,000  [45s]
done block a=652, up to 6,540,000,000  [46s]
done block a=653, up to 6,550,000,000  [46s]
done block a=654, up to 6,560,000,000  [46s]
done block a=655, up to 6,570,000,000  [46s]
done block a=656, up to 6,580,000,000  [46s]
done block a=657, up to 6,590,000,000  [47s]
done block a=658, up to 6,600,000,000  [47s]
done block a=659, up to 6,610,000,000  [47s]
done block a=660, up to 6,620,000,000  [47s]
done block a=661, up to 6,630,000,000  [48s]
done block a=662, up to 6,640,000,000  [48s]
done block a=663, up to 6,650,000,000  [48s]
done block a=664, up to 6,660,000,000  [48s]
done block a=665, up to 6,670,000,000  [49s]
done block a=666, up to 6,680,000,000  [49s]
done block a=667, up to 6,690,000,000  [49s]
done block a=668, up to 6,700,000,000  [49s]
done block a=669, up to 6,710,000,000  [50s]
done block a=670, up to 6,720,000,000  [50s]
done block a=671, up to 6,730,000,000  [50s]
done block a=672, up to 6,740,000,000  [50s]
done block a=673, up to 6,750,000,000  [51s]
done block a=674, up to 6,760,000,000  [51s]
done block a=675, up to 6,770,000,000  [51s]
done block a=676, up to 6,780,000,000  [51s]
done block a=677, up to 6,790,000,000  [52s]
done block a=678, up to 6,800,000,000  [52s]
done block a=679, up to 6,810,000,000  [52s]
done block a=680, up to 6,820,000,000  [52s]
done block a=681, up to 6,830,000,000  [53s]
done block a=682, up to 6,840,000,000  [53s]
done block a=683, up to 6,850,000,000  [53s]
done block a=684, up to 6,860,000,000  [53s]
done block a=685, up to 6,870,000,000  [53s]
done block a=686, up to 6,880,000,000  [53s]
done block a=687, up to 6,890,000,000  [54s]
done block a=688, up to 6,900,000,000  [54s]
done block a=689, up to 6,910,000,000  [54s]
done block a=690, up to 6,920,000,000  [54s]
done block a=691, up to 6,930,000,000  [55s]
done block a=692, up to 6,940,000,000  [55s]
done block a=693, up to 6,950,000,000  [55s]
done block a=694, up to 6,960,000,000  [55s]
done block a=695, up to 6,970,000,000  [56s]
done block a=696, up to 6,980,000,000  [56s]
done block a=697, up to 6,990,000,000  [56s]
done block a=698, up to 7,000,000,000  [56s]
done block a=699, up to 7,010,000,000  [57s]
done block a=700, up to 7,020,000,000  [57s]
done block a=701, up to 7,030,000,000  [57s]
done block a=702, up to 7,040,000,000  [57s]
done block a=703, up to 7,050,000,000  [57s]
done block a=704, up to 7,060,000,000  [58s]
done block a=705, up to 7,070,000,000  [58s]
done block a=706, up to 7,080,000,000  [58s]
done block a=707, up to 7,090,000,000  [58s]
done block a=708, up to 7,100,000,000  [59s]
done block a=709, up to 7,110,000,000  [59s]
done block a=710, up to 7,120,000,000  [59s]
done block a=711, up to 7,130,000,000  [59s]
done block a=712, up to 7,140,000,000  [59s]
done block a=713, up to 7,150,000,000  [60s]
done block a=714, up to 7,160,000,000  [60s]
done block a=715, up to 7,170,000,000  [60s]
done block a=716, up to 7,180,000,000  [60s]
done block a=717, up to 7,190,000,000  [61s]
done block a=718, up to 7,200,000,000  [61s]
done block a=719, up to 7,210,000,000  [61s]
done block a=720, up to 7,220,000,000  [61s]
done block a=721, up to 7,230,000,000  [62s]
done block a=722, up to 7,240,000,000  [62s]
done block a=723, up to 7,250,000,000  [62s]
done block a=724, up to 7,260,000,000  [62s]
done block a=725, up to 7,270,000,000  [63s]
done block a=726, up to 7,280,000,000  [63s]
done block a=727, up to 7,290,000,000  [63s]
done block a=728, up to 7,300,000,000  [63s]
done block a=729, up to 7,310,000,000  [63s]
done block a=730, up to 7,320,000,000  [64s]
done block a=731, up to 7,330,000,000  [64s]
done block a=732, up to 7,340,000,000  [64s]
done block a=733, up to 7,350,000,000  [64s]
done block a=734, up to 7,360,000,000  [65s]
done block a=735, up to 7,370,000,000  [65s]
done block a=736, up to 7,380,000,000  [65s]
done block a=737, up to 7,390,000,000  [65s]
done block a=738, up to 7,400,000,000  [66s]
done block a=739, up to 7,410,000,000  [66s]
done block a=740, up to 7,420,000,000  [66s]
done block a=741, up to 7,430,000,000  [66s]
done block a=742, up to 7,440,000,000  [66s]
done block a=743, up to 7,450,000,000  [67s]
done block a=744, up to 7,460,000,000  [67s]
done block a=745, up to 7,470,000,000  [67s]
done block a=746, up to 7,480,000,000  [67s]
done block a=747, up to 7,490,000,000  [68s]
done block a=748, up to 7,500,000,000  [68s]
done block a=749, up to 7,510,000,000  [68s]
done block a=750, up to 7,520,000,000  [68s]
done block a=751, up to 7,530,000,000  [68s]
done block a=752, up to 7,540,000,000  [69s]
done block a=753, up to 7,550,000,000  [69s]
done block a=754, up to 7,560,000,000  [69s]
done block a=755, up to 7,570,000,000  [70s]
done block a=756, up to 7,580,000,000  [70s]
done block a=757, up to 7,590,000,000  [70s]
done block a=758, up to 7,600,000,000  [70s]
done block a=759, up to 7,610,000,000  [70s]
done block a=760, up to 7,620,000,000  [71s]
done block a=761, up to 7,630,000,000  [71s]
done block a=762, up to 7,640,000,000  [71s]
done block a=763, up to 7,650,000,000  [71s]
done block a=764, up to 7,660,000,000  [71s]
done block a=765, up to 7,670,000,000  [72s]
done block a=766, up to 7,680,000,000  [72s]
done block a=767, up to 7,690,000,000  [72s]
done block a=768, up to 7,700,000,000  [72s]
done block a=769, up to 7,710,000,000  [73s]
done block a=770, up to 7,720,000,000  [73s]
done block a=771, up to 7,730,000,000  [73s]
done block a=772, up to 7,740,000,000  [73s]
done block a=773, up to 7,750,000,000  [74s]
done block a=774, up to 7,760,000,000  [74s]
done block a=775, up to 7,770,000,000  [74s]
done block a=776, up to 7,780,000,000  [74s]
done block a=777, up to 7,790,000,000  [74s]
done block a=778, up to 7,800,000,000  [75s]
done block a=779, up to 7,810,000,000  [75s]
done block a=780, up to 7,820,000,000  [75s]
done block a=781, up to 7,830,000,000  [75s]
done block a=782, up to 7,840,000,000  [76s]
done block a=783, up to 7,850,000,000  [76s]
done block a=784, up to 7,860,000,000  [76s]
done block a=785, up to 7,870,000,000  [76s]
done block a=786, up to 7,880,000,000  [76s]
done block a=787, up to 7,890,000,000  [77s]
done block a=788, up to 7,900,000,000  [77s]
done block a=789, up to 7,910,000,000  [77s]
done block a=790, up to 7,920,000,000  [78s]
done block a=791, up to 7,930,000,000  [78s]
done block a=792, up to 7,940,000,000  [78s]
done block a=793, up to 7,950,000,000  [78s]
done block a=794, up to 7,960,000,000  [79s]
done block a=795, up to 7,970,000,000  [79s]
done block a=796, up to 7,980,000,000  [79s]
done block a=797, up to 7,990,000,000  [79s]
done block a=798, up to 8,000,000,000  [79s]
Checked up to 8,000,000,000 in 79.4s
No exceptions found.
"""
