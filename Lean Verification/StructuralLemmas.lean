import Mathlib
import Defs

/-!
# Structural Lemmas

This file formalizes the key structural reductions used in the proof of Theorem 5.1,
which do not depend on analytic number theory or computation.

## Main results

- `hasRepresentation_iff_Tk_pos`: A number has a representation iff Tk > 0
- `excluded_set_properties`: Properties of the excluded k values
- `even_k_reduction`: Reduction from even k to odd m
- `representation_from_Tk_minus_Bq`: Combinatorial core of the analytic lemmas
-/

open Finset Nat

/-! ### Basic properties of representations -/

/-- If p is prime, p < n, n - p is squarefree and coprime to k,
    then n has a representation. -/
theorem hasRepresentation_of_prime_squarefree_coprime
    {k n p : ℕ} (hp : p.Prime) (hpn : p < n)
    (hsq : Squarefree (n - p)) (hcop : Nat.Coprime (n - p) k) :
    HasRepresentation k n :=
  ⟨p, hp, hpn, hsq, hcop⟩

/-- If k divides k', any representation coprime to k' is also coprime to k. -/
theorem hasRepresentation_of_dvd {k k' n : ℕ} (hd : k ∣ k')
    (h : HasRepresentation k' n) : HasRepresentation k n := by
  obtain ⟨p, hp, hpn, hsq, hcop⟩ := h
  exact ⟨p, hp, hpn, hsq, hcop.coprime_dvd_right hd⟩

/-- Representation with k = 1 is just prime + squarefree (Dudek's result). -/
theorem hasRepresentation_one_iff {n : ℕ} :
    HasRepresentation 1 n ↔ ∃ p, p.Prime ∧ p < n ∧ Squarefree (n - p) := by
  constructor
  · rintro ⟨p, hp, hpn, hsq, _⟩
    exact ⟨p, hp, hpn, hsq⟩
  · rintro ⟨p, hp, hpn, hsq⟩
    exact ⟨p, hp, hpn, hsq, Nat.coprime_one_right _⟩

/-! ### Properties of the excluded set -/

/-- 105 = 3 × 5 × 7 -/
theorem excluded_105 : 105 = 3 * 5 * 7 := by norm_num

/-- 429 = 3 × 11 × 13 -/
theorem excluded_429 : 429 = 3 * 11 * 13 := by norm_num

/-- All elements of the excluded set are squarefree -/
theorem excluded_squarefree : ∀ k ∈ ExcludedSet, Squarefree k := by
  intro k hk
  simp [ExcludedSet] at hk
  rcases hk with rfl | rfl <;> native_decide

/-- All elements of the excluded set are odd -/
theorem excluded_odd : ∀ k ∈ ExcludedSet, ¬ 2 ∣ k := by
  intro k hk
  simp [ExcludedSet] at hk
  rcases hk with rfl | rfl <;> omega

/-! ### Even k reduction (proof of Theorem 5.1, even case) -/

/-- For even n and odd prime p ≠ 2, n - p is odd and hence not divisible by 2. -/
theorem even_sub_odd_prime_odd {n p : ℕ} (hn : 2 ∣ n) (hp : Nat.Prime p)
    (hp2 : p ≠ 2) (hpn : p < n) : ¬ 2 ∣ (n - p) := by
  cases Nat.Prime.eq_two_or_odd hp <;> omega

/-- A number `n` has a representation as odd prime + squarefree coprime to `k`. -/
def HasOddRepresentation (k n : ℕ) : Prop :=
  ∃ p, p.Prime ∧ p ≠ 2 ∧ p < n ∧ Squarefree (n - p) ∧ Nat.Coprime (n - p) k

/-- Reduction from even k to odd m = k/2 for even n.
    For even k = 2m, if n is even and n = p + η with p an *odd* prime,
    η squarefree, (η, m) = 1, then η is odd so (η, 2) = 1 hence (η, k) = 1.

    Note: We need the representation to use an odd prime, since if p = 2
    and n is even, then n - 2 is even so not coprime to 2m. -/
theorem even_k_reduction {k m n : ℕ} (hk : k = 2 * m) (_hm_odd : ¬ 2 ∣ m)
    (hn_even : 2 ∣ n) (h : HasOddRepresentation m n) :
    HasRepresentation k n := by
  obtain ⟨p, hp, hp2, hpn, hsq, hcop⟩ := h
  exact ⟨p, hp, hpn, hsq, by
    simpa [hk] using Nat.Coprime.mul_right
      (show Nat.Coprime (n - p) 2 from
        Nat.Coprime.symm <| Nat.prime_two.coprime_iff_not_dvd.mpr <|
          even_sub_odd_prime_odd hn_even hp hp2 hpn)
      hcop⟩

/-! ### Key structural identity: R_k(n) and B_q(n) -/

/-- `Bq q n` counts primes p < n such that q ∣ (n - p) and n - p is squarefree.
    This is a simplified version of B_q(n) from the proof of Lemma 5.9. -/
def Bq (q n : ℕ) : ℕ :=
  ((Finset.range n).filter (fun p =>
    p.Prime ∧ q ∣ (n - p) ∧ Squarefree (n - p) ∧ p < n)).card

/-- The key structural bound: if T_{q₁}(n) > B_{q₂}(n) + B_{q₃}(n),
    then there exists a representation coprime to k = q₁ * q₂ * q₃.
    This is the combinatorial core of Lemmas 5.9–5.11. -/
theorem representation_from_Tk_minus_Bq
    {q₁ q₂ q₃ n : ℕ}
    (hq₁ : q₁.Prime) (hq₂ : q₂.Prime) (hq₃ : q₃.Prime)
    (h12 : q₁ ≠ q₂) (h13 : q₁ ≠ q₃) (_h23 : q₂ ≠ q₃)
    (_hn1 : Nat.Coprime n q₁) (_hn2 : Nat.Coprime n q₂) (_hn3 : Nat.Coprime n q₃)
    (hbound : Bq q₂ n + Bq q₃ n < Tk q₁ n) :
    HasRepresentation (q₁ * q₂ * q₃) n := by
  have h_exists_rep : ∃ p, p.Prime ∧ p < n ∧ Squarefree (n - p) ∧
      Nat.Coprime (n - p) q₁ ∧ ¬q₂ ∣ (n - p) ∧ ¬q₃ ∣ (n - p) := by
    contrapose! hbound
    refine le_trans (Finset.card_le_card ?_) (Finset.card_union_le _ _)
    grind +ring
  obtain ⟨p, hp₁, hp₂, hp₃, hp₄, hp₅, hp₆⟩ := h_exists_rep
  use p
  simp_all +decide [Nat.coprime_mul_iff_right]
  exact ⟨⟨hp₄, Nat.Coprime.symm <| hq₂.coprime_iff_not_dvd.mpr hp₅⟩,
    Nat.Coprime.symm <| hq₃.coprime_iff_not_dvd.mpr hp₆⟩
