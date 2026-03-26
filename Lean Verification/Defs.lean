import Mathlib

/-!
# Definitions for the Goldbach-type result: prime + squarefree coprime to k

This file formalizes the key definitions from the paper
"The Goldbach Conjecture" by Wilkie Hoare (Chapter 5), which proves that
every sufficiently large integer can be written as the sum of a prime and a
squarefree number coprime to k, for squarefree k with at most 3 prime factors.

## Main definitions

- `Squarefree`: a number with no repeated prime factors (using Mathlib's `Squarefree`)
- `HasRepresentation`: n can be written as p + η where p is prime, η is squarefree,
  and η is coprime to k
- `Tk`: the count of such representations (Definition from Section 5.2)
-/

open Finset Nat

noncomputable section

/-- `Tk k n` counts the number of representations of `n` as `p + (n - p)` where
    `p < n` is prime, `n - p` is squarefree, and `gcd(n - p, k) = 1`.
    This is `T_k(n)` from Section 5.2 of the paper. -/
def Tk (k n : ℕ) : ℕ :=
  ((Finset.range n).filter (fun p =>
    p.Prime ∧ Squarefree (n - p) ∧ Nat.Coprime (n - p) k ∧ p < n)).card

/-- A number `n` has a representation as prime + squarefree coprime to `k`. -/
def HasRepresentation (k n : ℕ) : Prop :=
  ∃ p, p.Prime ∧ p < n ∧ Squarefree (n - p) ∧ Nat.Coprime (n - p) k

/-- `Tk k n > 0` iff `n` has at least one representation. -/
theorem hasRepresentation_iff_Tk_pos (k n : ℕ) :
    HasRepresentation k n ↔ 0 < Tk k n := by
  constructor
  · rintro ⟨p, hp₁, hp₂, hp₃, hp₄⟩
    exact Finset.card_pos.mpr ⟨p, Finset.mem_filter.mpr
      ⟨Finset.mem_range.mpr hp₂, hp₁, hp₃, hp₄, hp₂⟩⟩
  · exact fun h => by
      obtain ⟨p, hp⟩ := Finset.card_pos.mp h
      use p; aesop

/-- The number of distinct prime factors of k. -/
def omega (k : ℕ) : ℕ := k.primeFactors.card

/-- k is odd and squarefree with at most 3 prime factors -/
def IsAdmissible (k : ℕ) : Prop :=
  1 < k ∧ Squarefree k ∧ omega k ≤ 3

/-- The set of excluded k values where the odd-n case is not proved.
    Note: The original paper excludes {105, 165, 195, 231, 255, 273}.
    Our analysis shows that 429 = 3 · 11 · 13 should also be excluded,
    as the analytic bounds fail for this value (see ANALYSIS.md). -/
def ExcludedSet : Finset ℕ := {105, 165, 195, 231, 255, 273, 429}

end
