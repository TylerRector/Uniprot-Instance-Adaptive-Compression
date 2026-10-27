/-!
# Feature-guarded correction theorem (Lean 4)

An indexed protein is a vector of residue symbols (represented as natural
numbers). If every mismatch at a protected coordinate is corrected, then
decode equals the target at every protected coordinate. The proof does *not*
assume biological functionality, and does not verify Python byte serialization.

`maskOr` formalizes mandatory site protection plus arbitrary optional edits.
-/

namespace FunctionGuard

abbrev Protein (n : Nat) := Fin n → Nat
abbrev Mask (n : Nat) := Fin n → Bool
abbrev Patches (n : Nat) := Fin n → Option Nat

def encodeSelected {n : Nat} (target reference : Protein n) (selected : Mask n) : Patches n :=
  fun i => if selected i then
    if target i = reference i then none else some (target i)
  else none

def decode {n : Nat} (reference : Protein n) (patches : Patches n) : Protein n :=
  fun i => match patches i with
    | none => reference i
    | some aa => aa

def maskOr {n : Nat} (guard optional : Mask n) : Mask n :=
  fun i => guard i || optional i

theorem selected_exact {n : Nat}
    (target reference : Protein n) (selected : Mask n) (i : Fin n)
    (hselected : selected i = true) :
    decode reference (encodeSelected target reference selected) i = target i := by
  by_cases hsame : target i = reference i
  · simp [decode, encodeSelected, hselected, hsame]
  · simp [decode, encodeSelected, hselected, hsame]

theorem protected_exact_with_optional {n : Nat}
    (target reference : Protein n) (guard optional : Mask n) (i : Fin n)
    (hprotected : guard i = true) :
    decode reference (encodeSelected target reference (maskOr guard optional)) i = target i := by
  apply selected_exact
  simp [maskOr, hprotected]

theorem unselected_is_reference {n : Nat}
    (target reference : Protein n) (selected : Mask n) (i : Fin n)
    (hnot : selected i = false) :
    decode reference (encodeSelected target reference selected) i = reference i := by
  simp [decode, encodeSelected, hnot]

theorem no_correction_for_match {n : Nat}
    (target reference : Protein n) (selected : Mask n) (i : Fin n)
    (hsame : target i = reference i) :
    encodeSelected target reference selected i = none := by
  by_cases hselected : selected i = true
  · simp [encodeSelected, hselected, hsame]
  · have hfalse : selected i = false := Bool.eq_false_iff.mpr hselected
    simp [encodeSelected, hfalse]

theorem selected_mismatch_has_correction {n : Nat}
    (target reference : Protein n) (selected : Mask n) (i : Fin n)
    (hselected : selected i = true) (hneq : target i ≠ reference i) :
    encodeSelected target reference selected i = some (target i) := by
  simp [encodeSelected, hselected, hneq]

theorem mismatch_requires_correction {n : Nat}
    (target reference : Protein n) (patches : Patches n) (i : Fin n)
    (hneq : reference i ≠ target i)
    (hexact : decode reference patches i = target i) :
    patches i ≠ none := by
  intro hnone
  have href : decode reference patches i = reference i := by
    simp [decode, hnone]
  exact hneq (href.symm.trans hexact)

#print axioms protected_exact_with_optional
#print axioms mismatch_requires_correction

end FunctionGuard
