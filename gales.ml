(*TO DO

UPDATES:
reconstruction power function to operate over zarith [o]

PROBABILITY COMPONENTS:
continuous distributions as probabilities, expectations, variances for
Gamma [o]
Beta [o]
Chi2 [o]

UTIL COMPONENTS:
graph module for random walks on graphs [o]
linear algebra components [o]

STOCHASTIC COMPONENTS:
RW on Graphs [o]
RW on Z2 (finish) [o]
GW processes [o]
Poisson processes [o]
Brownian motion [o]
MCMC [o]
Martingales [o]

STATISTICAL COMPONENTS:
TBD

*)

open List
open Float

(*LIST UTIL FUNCTIONS*)

(*applies a given function to two lists to return a chosen concatonation*)
let rec zipwith f l1 l2 =
  let rec zip_help f l1 l2 acc =
    match (l1,l2) with
    | (_,[]) -> acc
    | ([],_) -> acc
    | (h1::t1, h2::t2) -> (zip_help f t1 t2 (acc @ [f h1 h2]))
    in zip_help f l1 l2 []

(*returns list of numbers within a range*)
let rec range num1 num2 =
  if num2 < num1 then [ ]
  else num1::(range (num1 + 1) num2)

(*checks if a number exists in a list*)
let exists k l =
  List.fold_left (fun b x -> b || x = k) false l


(*MATH UTIL FUNCTIONS*)

(*recursive factorial function*)
let rec ffact (n : Z.t) (k : Z.t)=
  if Z.equal n k then Z.one
  else Z.mul n (ffact (Z.sub n Z.one) k)
  
(*cast a integer into Z*)
let cast n = Z.of_int n

(*recursively compute binomial coefficients*)
let binom_coef n k =
  let lg_n = cast n in let lg_k = cast k in
  let lg_d = Z.sub lg_n lg_k in
  (Q.make (ffact lg_n lg_d) (ffact lg_k Z.zero))
  
(*hard code "e" constant approximation instead of recursive estimation*)
let e_cnst =
  let e_num = (cast 2718281828459045235) in
  let e_den = (cast 1000000000000000000) in
  Q.make e_num e_den
  
(*hard code "pi" constant approximation instead of recursive estimation*)
let pi_cnst =
  let pi_num = (cast 3141592653589793238) in
  let pi_den = (cast 1000000000000000000) in
  Q.make pi_num pi_den

(*generation and handling of a sequence of IID random variables*)
(*let create_sq (x : RV.rv) (n : int) =
  let rec cnstrct_sq x n acc =
  match n with 
  | 0 -> Sq (acc)
  | _ -> cnstrct_sq x (n-1) (acc @ [x])
  in cnstrct_sq x n []*)

(*gamma function approximation using Lanczos method*)
let gamma (x : float) =
  if x < 1. then invalid_arg "Out of Range"
  else if x = 1. then Q.one
  else if ((Q.den (Q.of_float x)) = Z.one) then Q.make (ffact (cast (int_of_float (x-.1.))) Z.one) Z.one
  else
  let a = 10 in
  let z = x -. 1. in
  let rec gamma_build_sums k acc =
    match k with
    | 0 -> acc
    | _ -> 
      let sums_term1 = Q.inv (Q.of_float (z +. float_of_int k)) in
      let sums_term2 = Q.make (Z.( ** ) Z.minus_one (k-1)) (ffact (cast (k-1)) Z.zero) in
      let sums_term3 = Q.of_float (Float.pow (float_of_int (a - k)) ((float_of_int k) -. 0.5)) in
      let sums_term4 = Q.make (Z.( ** ) e_cnst.num (a - k)) (Z.( ** ) e_cnst.den (a - k)) in
      let sums_term_final = Q.( * ) sums_term1 (Q.( * ) sums_term2 (Q.( * ) sums_term3 sums_term4)) in
      gamma_build_sums (k - 1) (Q.(+) acc sums_term_final)
  in
  let sums = gamma_build_sums (a-1) (Q.of_float (Float.sqrt (2. *. (Q.to_float pi_cnst)))) in
  let term1 = Q.of_float (Float.pow (z +. (float_of_int a)) (z +. 0.5)) in
  let term2 = Q.of_float  (1. /. (Float.pow (Q.to_float e_cnst) (z +. (float_of_int a)))) in
  Q.( * ) term1 (Q.( * ) term2 sums)

(*use RV."name" to access components of module*)
module RV = struct
  type rv =
  | Disc_Unif of int
  | Bern of Q.t
  | RWBern of Q.t
  | Binom of int * Q.t
  | Geo of Q.t
  | Pois of Q.t
  | Cts_Unif of Q.t * Q.t
  | Exp of Q.t
  | Gauss of Q.t * Q.t
  | RWStopTime of rv
  | Sum of rv * rv
  
  (*test if valid probability*)
  let prob_meas p = if p <= 1.0 && p >= 0.0 then true else false

  (*create RVs by typing "let x = RV.'name' [necessary builders from below]"*)
  let disc_unif (n : int) =
    if n > 0 then Disc_Unif (n)
    else invalid_arg "Out of Range"
  
  let bern (p : float) =
    if prob_meas p then Bern (Q.of_float p)
    else invalid_arg "Out of Range"

  let rwbern (p : float) =
    if prob_meas p then RWBern (Q.of_float p)
    else invalid_arg "Out of Range"
  
  let binom (n : int) (p : float) =
    if n > 0 && prob_meas p then Binom (n, Q.of_float p)
    else invalid_arg "Out of Range"

  let geo (p : float) =
    if prob_meas p then Geo (Q.of_float p)
    else invalid_arg "Out of Range"
    
  let pois (lambda : float) =
    if lambda > 0.0 then Pois (Q.of_float lambda)
    else invalid_arg "Out of Range"

  let cts_unif (a : float) (b : float) =
    if a < b then Cts_Unif (Q.of_float a, Q.of_float b)
    else invalid_arg "Out of Range"

  let exp (lambda : float) =
    if lambda > 0.0 then Exp (Q.of_float lambda)
    else invalid_arg "Out of Range"
  
  let gauss (mu : float) (sig2 : float) =
    if sig2 >= 0.0 then Gauss (Q.of_float mu, Q.of_float sig2)
    else invalid_arg "Out of Range"

  let sum (x : rv) (n : int) =
    if n < 2 then invalid_arg "Out of Range"
    else let rec cnstrct_sum x k acc =
    match k with
    | 0 -> acc
    | _ -> cnstrct_sum x (k-1) (Sum (x, acc))
    in cnstrct_sum x (n-2) (Sum (x, x))

end

(*output discrete probably as float*)
let rec disc_pr (x : RV.rv) (value : float) = 
  let v = Q.of_float value in
  match x with
  
  | RV.Disc_Unif (n) ->
    if Z.equal v.den Z.one && exists (Q.to_int v) (range 1 n)
    then Q.to_float (Q.make Z.one (cast n))
    else 0.0
  
  | RV.Bern (p) -> 
    if Q.equal v Q.one then Q.to_float p
    else if Q.equal v Q.zero then Q.to_float (Q.(-) Q.one p)
    else 0.0
  
  | RV.RWBern (p) ->
    if Q.equal v Q.one then Q.to_float p
    else if Q.equal v Q.zero then Q.to_float (Q.(-) Q.one p)
    else 0.0
  
  | RV.Geo (p) ->
    let k = Z.to_int v.num in
    if k >= 0 && Z.equal v.den Z.one then
    let comp = Q.(-) (Q.make Z.one Z.one) p in
    let p2 = Q.make (Z.( ** ) comp.num k) (Z.( ** ) comp.den k) in
    Q.to_float (Q.( * ) p p2)
    else 0.0

  | RV.Binom (n, p) -> 
    let k = Z.to_int v.num in
    if k >= 0 && Z.equal v.den Z.one then
    let p1 = Q.make (Z.( ** ) p.num k) (Z.( ** ) p.den k) in
    let comp = Q.(-) (Q.make Z.one Z.one) p in
    let p2 = Q.make (Z.( ** ) comp.num (n-k)) (Z.( ** ) comp.den (n-k)) in
    Q.to_float (Q.( * ) (binom_coef n k) (Q.( * ) p1 p2))
    else 0.0

  | RV.Pois (lambda) ->
    let k = Z.to_int v.num in
    if k >= 0 && Z.equal v.den Z.one then
    let inv_e_lam = Q.inv (Q.of_float (Float.pow (Q.to_float e_cnst) (Q.to_float lambda))) in
    let lam_k = Q.make (Z.( ** ) (lambda.num) k) (Z.( ** ) (lambda.den) k) in
    Q.to_float (Q.mul (Q.mul lam_k inv_e_lam) (Q.make Z.one (ffact (cast k) Z.zero)))
    else 0.0

  | RV.Cts_Unif (a, b) -> invalid_arg "Not Discrete"
  
  | RV.Gauss (mu, sig2) -> invalid_arg "Not Discrete"

  | RV.Exp (lambda) -> invalid_arg "Not Discrete"

  | _ -> 0.0

(*output continuous probability as float*)
let rec cts_pr (x : RV.rv) (lower : float) (upper : float) =
  if lower > upper then invalid_arg "Out of Range" else
  let lb = Q.of_float lower in
  let ub = Q.of_float upper in
  match x with

  | RV.Disc_Unif (n) -> invalid_arg "Not Continuous"

  | RV.Bern (p) -> invalid_arg "Not Continuous"

  | RV.RWBern (p) -> invalid_arg "Not Continuous"

  | RV.Binom (n, p) -> invalid_arg "Not Continuous"

  | RV.Pois (lambda) -> invalid_arg "Not Continuous"

  | RV.Cts_Unif (a, b) ->
    let den_cnst = (Q.inv (Q.(-) b a)) in
    if Q.leq lb a && Q.leq b ub then 1.0
    else if Q.leq lb a && Q.leq ub b then
    Q.to_float (Q.( * ) den_cnst (Q.(-) ub a))
    else if Q.leq a lb && Q.leq b ub then
    Q.to_float (Q.( * ) den_cnst (Q.(-) b lb))
    else Q.to_float (Q.( * ) den_cnst (Q.(-) ub lb))

  | RV.Exp (lambda) ->
    if lower >= 0. && upper >= 0. then
    let e_lam_lb = Float.pow (Q.to_float e_cnst) (Q.to_float (Q.( * ) lambda lb)) in
    let e_lam_ub = Float.pow (Q.to_float e_cnst) (Q.to_float (Q.( * ) lambda ub)) in
    Q.to_float (Q.(-) (Q.inv (Q.of_float e_lam_lb)) (Q.inv (Q.of_float e_lam_ub)))
    else if lower < 0. && upper >= 0. then
    let e_lam_lb = Float.pow (Q.to_float e_cnst) (Q.to_float (Q.( * ) lambda Q.zero)) in
    let e_lam_ub = Float.pow (Q.to_float e_cnst) (Q.to_float (Q.( * ) lambda ub)) in
    Q.to_float (Q.(-) (Q.inv (Q.of_float e_lam_lb)) (Q.inv (Q.of_float e_lam_ub)))
    else 0.0

  | RV.Gauss (mu, sig2) ->
    let den_cnst = Q.inv (Q.of_float (Float.sqrt(Q.to_float (Q.( * ) (Q.make (cast 2) Z.one) sig2)))) in
    let ub_factor = Q.( * ) (Q.(-) ub mu) den_cnst in
    let lb_factor = Q.( * ) (Q.(-) lb mu) den_cnst in
    let diff = Q.(-) (Q.of_float (Float.erf (Q.to_float ub_factor))) (Q.of_float (Float.erf (Q.to_float lb_factor))) in 
    Q.to_float (Q.( * ) (Q.make Z.one (cast 2)) diff)

  | _ -> 0.0

(*output expectation as float*)
let rec e (x : RV.rv) =
  match x with
  
  | RV.Disc_Unif (n) ->
      let vals = List.map float_of_int (range 1 n) in
      let probs = (List.map (disc_pr x) (List.map float_of_int (range 1 n))) in
      List.fold_left (+.) 0.0 (zipwith ( *. ) vals probs)

  | RV.Bern (p) -> Q.to_float p

  | RV.RWBern (p) -> Q.to_float (Q.(-) (Q.( * ) (Q.make (cast 2) Z.one) p) Q.one)
  
  | RV.Binom (n, p) -> Q.to_float (Q.( * ) (Q.make (cast n) Z.one) p)

  | RV.Geo (p) -> Q.to_float (Q.( * ) (Q.(-) Q.one p) (Q.inv p))

  | RV.Pois (lambda) -> Q.to_float (lambda)

  | RV.Cts_Unif (a, b) -> Q.to_float (Q.( * ) (Q.add a b) (Q.make (cast 1) (cast 2)))

  | RV.Exp (lambda) -> Q.to_float (Q.inv lambda)

  | RV.Gauss (mu, sig2) -> Q.to_float (mu)

  | RV.Sum (z, tail) -> (e z) +. (e tail)

  | _ -> 0.0

(*output variance as float*)
let rec var (x : RV.rv) =
  match x with
  
  | RV.Disc_Unif (n) ->
    let n2 = Z.( ** ) (Z.of_int n) 2 in
    let diff = Q.(-) (Q.make n2 Z.one) Q.one in
    Q.to_float (Q.( * ) (Q.make Z.one (cast 12)) diff)

  | RV.Bern (p) -> Q.to_float (Q.( * ) p (Q.(-) Q.one p))
  
  | RV.Binom (n, p) -> Q.to_float (Q.( * ) (Q.( * ) (Q.make (cast n) Z.one) p) (Q.(-) Q.one p))

  | RV.Geo (p) ->
    let p2_num = Z.( ** ) p.num 2 in
    let p2_den = Z.( ** ) p.den 2 in
    let p2 = Q.make p2_num p2_den in
    Q.to_float (Q.( * ) (Q.(-) Q.one p) (Q.inv p2))

  | RV.Pois (lambda) -> Q.to_float (lambda)

  | RV.Cts_Unif (a, b) ->
    let diff = Q.(-) b a in
    let diff2_num = Z.( ** ) diff.num 2 in
    let diff2_den = Z.( ** ) diff.den 2 in
    let diff2 = Q.make diff2_num diff2_den in
    Q.to_float (Q.( * ) (Q.make Z.one (cast 12)) diff2)

  | RV.Exp (lambda) -> 
    let lambda2_num = Z.( ** ) lambda.num 2 in
    let lambda2_den = Z.( ** ) lambda.den 2 in
    let lambda2 = Q.make lambda2_num lambda2_den in
    Q.to_float (Q.inv lambda2)
    
  | RV.Gauss (mu, sig2) -> Q.to_float (sig2)

  | _ -> 0.0