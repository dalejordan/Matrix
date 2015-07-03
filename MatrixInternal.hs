-- Internal Matrix operations using mutable Vector operations for speed

{-# LANGUAGE BangPatterns, ScopedTypeVariables #-}
{-# OPTIONS -fno-warn-name-shadowing #-}

module MatrixInternal (luDecomp_,solve_) where

-- These are separate mainly to avoid import collisions

import Control.Monad
import Control.Monad.ST
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as M

-- LU Decomposition

luDecomp_ :: (Ord a,Fractional a) => V.Vector a -> Int -> Maybe (V.Vector a,V.Vector Int)
luDecomp_ v n = runST (lud v n)

lud :: (Ord a,Fractional a) => V.Vector a -> Int -> ST s (Maybe (V.Vector a,V.Vector Int))
lud v n = do
  p :: M.STVector s Int <- M.new (n-1)
  s :: M.STVector s Int <- makePV n
  x :: M.STVector s a <- M.new (V.length v)
  V.unsafeCopy x v
  ok <- loopM 0 (n-1) True $ \i ok -> 
         if ok
            then do
              i' <- (n*) <$> M.read s i
              xmx0 <- abs <$> M.read x (i'+i)
              (xmx,j) <- loopM (i+1) n (xmx0,0) $ \k (xmx,j) -> do
                xmx' <- abs <$> M.read x (k*n+i)
                return $! if xmx >= xmx' then (xmx,j) else (xmx',k)
              if xmx == 0
                 then return False
                 else do
                   M.write p i j
                   when (i /= j) $ do
                     M.swap s i j
                     let j' = n*j
                     loopM_ i n $ \k -> M.swap x (i'+k) (j'+k)
                   xii <- M.read x (i'+i)
                   loopM_ (i+1) n $ \r -> do
                     let r' = n*r
                     xri <- (/xii) <$> M.read x (r'+i)
                     M.write x (r'+i) $ xri
                     loopM_ (i+1) n $ \c -> do
                       xrc <- M.read x (r'+c)
                       xic <- M.read x (i'+c)
                       M.write x (r'+c) $! xrc - xri*xic
                   return True
            else return False
  xnn <- M.read x (n*n-1)
  if ok && xnn /= 0
     then do
       x' <- V.freeze x
       p' <- V.freeze p
       return $! Just (x',p')
     else return Nothing

-- General linear solver, called for various end uses

solve_ :: (Fractional a)
       => V.Vector a                    -- Constant value matrix (n*q) (column major order)
       -> (V.Vector a,V.Vector Int)      -- LU decomposition matrix (n*n) and pivot vector
       -> V.Vector a                     -- Result matrix (n*q) (column major order)
solve_ b (lu,p) = V.create $ do
  let n = V.length p + 1
      cols = V.length b `div` n
  x :: M.STVector s a <- M.new (V.length b)
  y :: M.STVector s a <- M.new n
  loopM_ 0 cols $ \c -> do
    let c' = c*n
    -- copy the column to y
    V.unsafeCopy y (V.slice c' n b)
    -- Forward substitution solving: L |*| pb = y
    loopM_ 0 (n-1) $ \i -> do
      let j = p V.! i
      when (i /= j) $ M.swap y i j
      yi <- M.read y i
      loopM_ (i+1) n $ \j -> do
        yj <- M.read y j
        M.write y j $! yj - yi * lu V.! (j*n+i)
    -- Backward substitution solving U |*| y = x
    yn <- M.read y (n-1)
    M.write x (cols*(n-1)+c) (yn / lu V.! (n*n-1))
    loopM_ 1 n $ \r ->
      let i = n - r - 1
          i' = i*cols
          i'' = i*n
          uii = lu V.! (i''+i)
      in do yi <- M.read y i
            z <- loopM (i+1) n 0 $ \j a -> do
                  let uij = lu V.! (i''+j)
                  xj <- M.read x (c+j*cols)
                  return $! a + uij * xj
            M.write x (i'+c) $! (yi - z)/uii
  return x

-- create permutation vector for pivot bookkeeping

makePV :: Int -> ST s (M.STVector s Int)
makePV n = do
  p <- M.new n
  loopM_ 0 (M.length p) $ \i -> M.unsafeWrite p i i
  return p

-- Simple monadic for loops

-- In ghc-7.8.2 and earlier, the usual idiom "forM_ [n0..n-1]" has
-- problems with not being optimized when n is a constant, due to the
-- list being floated out.

loopM_ :: (Monad m) => Int -> Int -> (Int -> m ()) -> m ()
{-# INLINE loopM_ #-}
loopM_ n0 n f = go n0
  where go !i | i == n     = return ()
              | otherwise = f i >> go (i+1)

loopM :: (Monad m) => Int -> Int -> a  -> (Int -> a -> m a) -> m a
{-# INLINE loopM #-}
loopM n0 n a0 f = go n0 a0
  where go !i !a | i == n   = return a
                 | otherwise = f i a >>= go (i+1)
