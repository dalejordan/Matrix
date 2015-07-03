-- Some simple matrix functions over Data.Vector

{-# OPTIONS -fno-warn-name-shadowing #-}

module Matrix
  (Matrix (M)
  ,shape
  ,theVector
  ,theRow
  ,theCol
  ,rowAt
  ,colAt
  ,(!)
  ,matrix
  ,fromList
  ,zero
  ,ones
  ,eye
  ,tp
  ,(|*|)
  ,inv
  ,inv'
  ,luDecomp
  ,luDecomp'
--  ,ppMatrix
--  ,prM
--  ,prIM
  ) where

import MatrixInternal
--import String

import Control.DeepSeq
import Data.Bifunctor
import Data.Maybe
import Data.Vector (Vector)
import qualified Data.Vector as V

data Matrix a = M
  { _nrows :: !Int                -- Number of rows.
  , _ncols :: !Int                -- Number of columns.
  , vec    :: Vector a            -- Content of the matrix
  } deriving (Eq,Show,Read)

instance Functor Matrix where
  {-# INLINE fmap #-}
  fmap f (M n m v) = M n m $ V.map f v

instance Foldable Matrix where
  foldMap f = foldMap f . vec

instance Traversable Matrix where
  sequenceA (M n m v) = fmap (M n m) $ sequenceA v

instance NFData a => NFData (Matrix a) where
  rnf (M _ _ v) = rnf v

instance Num a => Num (Matrix a) where
  fromInteger = M 1 1 . V.singleton . fromInteger
  negate = fmap negate
  abs = fmap abs
  signum = fmap signum
  (M n m v) + (M n' m' v')
    | n == n' && m == m' = M n m $ V.zipWith (+) v v'
    | otherwise       = error "Matrix.Add size mismatch"
  (M n m v) - (M n' m' v')
    | n == n' && m == m' = M n m $ V.zipWith (-) v v'
    | otherwise       = error "Matrix.Subtract size mismatch"
  (M n m v) * (M n' m' v')
    | n == n' && m == m' = M n m $ V.zipWith (*) v v'
    | otherwise       = error "Matrix.Multiply size mismatch"

instance Fractional a => Fractional (Matrix a) where
  (M n m v) / (M n' m' v')
    | n == n' && m == m' = M n m $ V.zipWith (/) v v'
    | otherwise       = error "Matrix.Divide size mismatch"
  recip = undefined              -- inv has an Ord constraint...
  fromRational = M 1 1 . V.singleton . fromRational

shape :: Matrix a -> (Int,Int)
shape (M n m _) = (n,m)

theVector :: Matrix a -> Vector a
theVector (M _ _ v) = v

theRow :: Matrix a -> Vector a
theRow (M 1 _ v) = v
theRow _         = error "theRow: not a row vector"

theCol :: Matrix a -> Vector a
theCol (M _ 1 v) = v
theCol _         = error "theCol: not a column vector"

rowAt :: Matrix a -> Int -> Vector a
rowAt (M m _ v) r = V.slice (encode m (r,1)) m v

colAt :: Matrix a -> Int -> Vector a
colAt (M m n v) c = V.generate n eleAt
  where eleAt i = v V.! encode m (i,c)

-- subscript handling

encode :: Int -> (Int,Int) -> Int
{-# INLINE encode #-}
encode m (i,j) = (i-1)*m + j - 1

decode :: Int -> Int -> (Int,Int)
{-# INLINE decode #-}
decode m k = (q+1,r+1)
  where (q,r) = quotRem k m

(!) :: Matrix a -> (Int,Int) -> a
{-# INLINE (!) #-}
(M _ m v) ! ij = v V.! encode m ij

-- utility constructors

matrix :: Int -> Int -> ((Int,Int) -> a) -> Matrix a
{-# INLINE matrix #-}
matrix n m f = M n m . V.generate (n*m) $ f . decode m

fromList :: Int -> Int -> [a] -> Matrix a
{-# INLINE fromList #-}
fromList n m = M n m . V.fromListN (n*m)

-- Some common "constant" matrices

zero :: Num a => Int -> Int -> Matrix a
zero n m = M n m $ V.replicate (n*m) 0

ones :: Num a => Int -> Int -> Matrix a
ones n m = M n m $ V.replicate (n*m) 1

eye :: Num a => Int -> Matrix a
eye n = matrix n n $ \(i,j) -> if i == j then 1 else 0

-- Transpose

tp :: Matrix a -> Matrix a
tp (M 1 m v) = M m 1 v
tp (M n 1 v) = M 1 n v
tp x@(M n m _) = matrix n m $ \(i,j) -> x ! (j,i)

-- Matrix Multiply

(|*|) :: Num a => Matrix a -> Matrix a -> Matrix a
{-# INLINE (|*|) #-}
(|*|) a1@(M n m _) a2@(M n' m' _)
    | m == n'    = matrix n m' $ \(i,j) -> sum [ a1 ! (i,k) * a2 ! (k,j) | k <- [1 .. m] ]
    | otherwise = error "Matrix.|*|: sizes invalid"

-- Matrix Inverse

-- Returns Nothing if matrix is singular

inv' :: (Ord a,Fractional a) => Matrix a -> Maybe (Matrix a)
inv' (M n m v) | n == m     = M n n . solve_ eye' <$> luDecomp_ v n
               | otherwise = error "Matrix.inv: non-sqaure matrix"
  where (M _ _ eye') = eye n

-- Errors if matrix is singular

inv :: (Ord a,Fractional a) => Matrix a -> Matrix a
inv = fromMaybe (error "Matrix.inv': singular matrix") . inv'

-- LU decomposition (with partial pivoting)

luDecomp' :: (Ord a,Fractional a) => Matrix a -> Maybe (Matrix a,Matrix Int)
luDecomp' (M n m v) | n == m     = bimap (M n n) (M (n-1) 1) <$> luDecomp_ v n
                    | otherwise = error "Matrix.luDecomp: non-square matrix"

luDecomp :: (Ord a,Fractional a) => Matrix a -> (Matrix a,Matrix Int)
luDecomp = fromMaybe (error "Matrix.luDecomp: singular matrix") . luDecomp'

-- Simple pretty printer
{-
ppMatrix :: (a -> String) -> Matrix a -> [String]
ppMatrix ppEle (M _ m v) = map ppVec $ slices v
  where slices v = map (\i -> V.slice i m v) starts
          where lv = V.length v
                starts = [0,m..lv-1]
        ppVec v = "| " ++ concatMap ((++ " ") . ppEle) (V.toList v) ++ "|"

prM :: (RealFloat a) => Matrix a -> IO ()
prM = mapM_ putStrLn . ppMatrix (justifyRight 7 . showFixed 3)

prIM :: Matrix Int -> IO ()
prIM = mapM_ putStrLn . ppMatrix (justifyRight 3 . show)
-}
