{-# LANGUAGE BangPatterns #-}
module Math.Sequence.Trig where

-- TODO: clean up this comment, I got a bit long-winded.

-- |Efficient recurrence for cosine and sine sequences with regularly-spaced
-- parameters.  @trigSeq theta delta@ === @zip (map sin thetas) (map cos thetas)@
-- where @thetas = iterate (+delta) theta@.
-- 
-- This recurrence is stable, but not always appropriate.  In typical situations,
-- error is dominated by drift in the original enumeration (for example,
-- in @[0,1/3 ..]@, 1/3 is inexact and the drift from this minor error will 
-- propagate differently in the enumeration than in this recurrence, leading 
-- to drift when comparing the expressions above for long sequences).  If the
-- original enumeration is exact (eg, @[0,0.125..]@) then the drift is quite 
-- slow, with a typical maximum error on the order of 1e-13 over a list of 
-- 100000 'Double's (compared to 1e-10 in the case where the enumeration is 
-- computed at higher precision, or 1e-7 when the enumeration is computed naively
-- - note that in the latter case the error is in the enumeration, not the
-- recurrence).
trigSeq :: Floating a => a -> a -> [(a,a)]
trigSeq theta delta = go (sin theta) (cos theta)
    where
        alpha = 2 * sin (delta / 2) ^ 2
        beta  = sin delta
        
        go !sin_theta !cos_theta
            = (sin_theta, cos_theta)
            : go next_sin next_cos
            where
                next_sin = sin_theta - (alpha * sin_theta - beta * cos_theta)
                next_cos = cos_theta - (alpha * cos_theta + beta * sin_theta)

-- |@enumSinFromBy theta delta@ === @map sin (iterate (+delta) theta)@
enumSinFromBy :: (Floating b) => b -> b -> [b]
enumSinFromBy theta delta = map fst (trigSeq theta delta)
-- |@enumCosFromBy theta delta@ === @map cos (iterate (+delta) theta)@
enumCosFromBy :: (Floating b) => b -> b -> [b]
enumCosFromBy theta delta = map snd (trigSeq theta delta)
-- |@enumTanFromBy theta delta@ === @map tan (iterate (+delta) theta)@
enumTanFromBy :: (Floating b) => b -> b -> [b]
enumTanFromBy theta delta = map (uncurry (/)) (trigSeq theta delta)