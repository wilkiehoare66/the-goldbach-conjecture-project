# The following code deals with odd n only, even n are verified using the is_good function under Lemma4.3.py.

from sympy import prevprime, primerange, nextprime, factorint, isprime

OMEGA_K_MAX = 35

def coprime_weaker(a,b):
    fi_A = factorint(a)
    fi_B = factorint(b)
    for p in fi_A:
        if p in fi_B.keys() and p > 2:
            return False
    return True

def nextcoprimesquarefree_evens(y,coprime_to):
    i = 2
    while i > 0:
        fi = factorint(y + i)
        condition_met = [coprime_weaker(j,y+i) for j in coprime_to] 
        if max(fi.values()) == 1 and False not in condition_met:
            return y + i
        i += 2

def is_good(n,v):
    sf = 2 #the first nontrivial squarefree
    cpt = []
    while sf < n:
        if isprime(n-sf) == True:
            cpt.append(sf)
        if len(cpt) < v+1:
            sf = nextcoprimesquarefree_evens(sf, cpt)
        else:
            return True, cpt
    return False, cpt

def is_prime_times_power_of_two(n):
    while n % 2 == 0:
        n = n/2
    return isprime(n)

def issquarefree(y):
    fi = factorint(y)
    if max(fi.values()) == 1:
        return True
    else:
        return False

def is_sqarefree_times_power_of_two(n):
    while n % 2 == 0:
        n = n/2
    return issquarefree(n)

def verify_ft_semiprimes_small(from_me, to_me, N=3):
    not_verified = []
    m = from_me
    while m <= to_me+1:
        ig = is_good(m,N)
        if ig[0] == False:
            #print(m)
            not_verified.append(m)
        m += 2
    return not_verified

def verify_ft_semiprimes_quick(from_me, to_me):
    small_primes = set()
    p = 2
    while p < 500:
        small_primes.add(p)
        p = nextprime(p)

    """
    pc_primes = set()
    p = prevprime(from_me - 3*10**3)
    while p < to_me:
        pc_primes.add(p)
        p = nextprime(p+10**4)
    """
    
    num_representations = {}
    q = prevprime(from_me - 10**3)
    while q < to_me:
        for p1 in small_primes:
            for p2 in small_primes:
                if p1 == p2:
                    m = q + 2*p1
                    if m not in num_representations:
                        num_representations[m] = 1
                    else:
                        num_representations[m] += 1
                else:
                    m = q + p1 + p2
                    for a1, a2 in [(p1,p1), (q,p1), (q,p2)]:
                        #print(a1,a2,[(p1,p1), (q,p1), (q,p2)])
                        if is_prime_times_power_of_two(a1 + a2) == True:
                            if m not in num_representations:
                                num_representations[m] = 1
                            else:
                                num_representations[m] += 1
        #print(q, len(num_representations))
        q = nextprime(q + 10**3)
    #print("Generated!")
    
    not_verified = []
    m = from_me
    while m <= to_me:
        if m not in num_representations or num_representations[m] < 2:
            #print(m)
            reps = []
            q = 3
            while 2*q < m and len(reps) < 2:
                if is_prime_times_power_of_two(m-q) == True:
                    #print(q,m-q)
                    reps.append((q,m-q))
                q = nextprime(q)
            #print(m, reps)
            if len(reps) < 2:
                not_verified.append(m)
        #else:
        #    print(m)
        m += 2
    return not_verified


#exptons = verify_ft_semiprimes_small(5, 10**5, 2)
#for e in exptons:
    #print(e)
#print("\n")

exptons = verify_ft_semiprimes_small(5, 10**5, OMEGA_K_MAX)
for e in exptons:
    print(e)

"""
Output:

5
7
9
11
13
15
17
19
21
23
25
27
29
31
33
35
37
39
41
43
45
47
49
51
53
55
57
59
61
63
65
67
69
71
73
75
77
79
81
83
85
87
89
91
93
95
97
99
101
103
105
107
109
111
113
115
117
119
121
123
125
127
129
131
133
135
137
139
141
143
145
147
149
151
153
155
157
159
161
163
165
167
169
171
173
175
177
179
181
183
185
187
189
191
193
195
197
199
201
203
205
207
209
211
213
215
217
219
221
223
225
227
229
231
233
235
237
239
241
243
245
247
249
251
253
255
257
259
261
263
265
267
269
271
273
275
277
279
281
283
285
287
289
291
293
295
297
299
301
303
305
307
309
311
313
315
317
319
321
323
325
327
329
331
333
335
337
339
341
343
345
347
349
351
353
355
357
359
361
363
365
367
369
371
373
375
377
379
381
383
385
387
389
391
393
395
397
399
401
403
405
407
409
411
413
415
417
419
421
423
425
427
429
431
433
435
437
439
441
443
445
447
449
451
453
455
457
459
461
463
465
467
469
471
473
475
477
479
481
483
485
487
489
491
493
495
497
499
501
503
505
507
509
511
513
515
517
519
521
523
527
529
531
533
535
537
539
541
543
545
547
549
551
553
555
557
559
561
563
565
567
569
571
573
575
577
579
581
583
587
589
591
593
595
597
599
601
603
605
607
611
613
617
619
621
623
625
627
629
631
633
635
637
639
641
643
647
649
651
653
655
657
659
661
663
665
667
669
671
673
677
679
683
685
687
689
691
695
697
699
701
703
707
709
711
713
715
717
719
721
723
725
727
729
731
733
737
739
743
745
747
749
751
755
757
761
763
767
769
773
775
779
781
785
787
791
793
797
799
803
805
809
811
815
817
821
823
827
829
833
835
839
841
843
845
847
851
853
857
859
863
865
869
871
875
877
881
883
887
889
893
895
899
901
905
907
911
913
917
919
923
925
929
931
935
937
941
943
947
949
953
955
959
961
965
967
971
973
977
979
983
985
989
991
995
997
1001
1003
1007
1009
1013
1019
1021
1025
1027
1031
1033
1037
1039
1043
1045
1049
1051
1055
1057
1061
1063
1067
1069
1073
1075
1079
1081
1087
1091
1093
1097
1099
1103
1105
1109
1111
1115
1117
1121
1123
1127
1129
1133
1135
1139
1141
1145
1147
1151
1153
1157
1159
1163
1165
1169
1171
1175
1177
1181
1183
1187
1189
1193
1195
1199
1201
1207
1211
1213
1217
1219
1223
1229
1231
1237
1241
1243
1247
1249
1253
1259
1261
1267
1271
1273
1277
1279
1283
1285
1289
1291
1297
1301
1303
1307
1313
1319
1321
1327
1331
1333
1337
1339
1343
1349
1351
1355
1357
1361
1363
1367
1369
1373
1379
1381
1387
1391
1393
1397
1399
1403
1409
1411
1417
1423
1427
1429
1433
1439
1441
1447
1451
1453
1457
1459
1465
1469
1471
1477
1481
1483
1487
1489
1493
1499
1501
1507
1511
1517
1523
1529
1531
1535
1537
1541
1543
1549
1559
1567
1571
1579
1583
1591
1597
1601
1607
1609
1613
1619
1621
1627
1637
1643
1651
1657
1663
1667
1669
1679
1691
1693
1697
1699
1703
1709
1711
1721
1727
1733
1739
1741
1753
1763
1783
1799
1801
1811
1819
1823
1831
1847
1867
1873
1879
1951
1979
1993
1999
2011
2027
2089
2183
"""
