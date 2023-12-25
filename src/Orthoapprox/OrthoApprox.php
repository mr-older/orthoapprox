<?php
/*
*	Polinomial regression via orthogonal regression polynomials
*	https://ru.stackoverflow.com/questions/383716/%D0%9C%D0%9D%D0%9A-%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC/639155
*	Author: Yury Negometyanov, 21.03.2017
*	Cut and edited for approximation needs by Evgeny Zakharenko, 24.12.2023
*
*	x and y are arrays of original dataset
*	valueByX returns approximated (by ortogonal polynomials) value
*
*/
class OrthoApprox
{
	private $x;
	private $y;
	private $m;
	private $values;
	private $sums;
	private $coeffs;

	public function __construct($x, $y, $m)
	{
		$this->x = $x;
		$this->y = $y;
		$this->m = $m;
		$this->calcValuesSums();
		$this->calcAllCoeffs();
	}

	public function subtract($ar, $subtr)
	{
		if(gettype($subtr) == "array") {
			return array_map(function($val1, $val2) { return $val1 - $val2; }, $ar, $subtr);
		}
		else return array_map(function($val) use($subtr) { return $val - $subtr; }, $ar);
	}

	public function multi($ar, $mult)
	{
		if($mult == "square") {
			return array_map(function($a) { return $a * $a; }, $ar);
		}
		elseif(gettype($mult) == "array") {
			return array_map(function($a, $b) { return $a * $b; }, $ar, $mult);
		}
		else return array_map(function($val) use ($mult) { return $val * $mult; }, $ar);
	}

	private function calcValuesSums()
	{
		for($j = -1; $j < $this->m; $j++) {
			if(!isset($v)) {
				$v = [array_fill(0, count($this->x), 1.0)];
				$s = [(float)count($this->x)];
				$q = [array_sum($this->x)];
			} else {
				$v[] = $this->subtract($this->x, end($q) / end($s)); 
				if($j > 0) {
					$v[$j + 1] = $this->subtract($this->multi(end($v), prev($v)), 
												$this->multi(prev($v), end($s) / prev($s)));
				}
				$sq = $this->multi(end($v), "square");
				$s[] = array_sum($sq);
				$q[] = array_sum($this->multi($this->x, $sq));
			}
		}
		$this->values = $v;
		$this->sums = ['∑P²' => $s, '∑xP²;' => $q];
	}

	private function calcOrthoCoeffs()
	{
		$prev = null;
		$ab = array_map(function($ss, $qq) use (&$prev)
		{
			if(is_null($prev)) $bb = 0;
			else $bb = $ss / $prev;
			$prev = $ss;
			return [$qq / $ss, $bb];
		}, reset($this->sums), next($this->sums));

		foreach($ab as $k => list($aa, $bb)) {
			if(!isset($с)) $с = [[1.0]];
			$с[] = array_fill(0, $k + 2, 1.0);
			for($j = 0; $j <= $k; $j++) { 
				$с[$k + 1][$j] = ($j == 0 ? 0.0 : $с[$k][$j-1])
								- $с[$k][$j] * $aa
								- (($j == $k) ? 0.0 : $с[$k-1][$j] * $bb);
			}
		}
		return $с;
	}

	private function calcAllCoeffs()
	{
		$s = reset($this->sums);
		$b = [];

		foreach($this->values as $key => $v) {
			$b[] = array_sum($this->multi($this->y, $v)) / $s[$key];
		}

		$c = $this->calcOrthoCoeffs();
		$m = count($b) - 1;
		$a = array_fill(0, $m + 1, 0.0);
		foreach($a as $j => &$aj) {
			for($k = $j; $k <= $m; $k++) {
				$aj+= $b[$k] * $c[$k][$j];
			}
		}

		$this->coeffs = ['c' => $c, 'b' => $b, 'a' => $a];
	}

	public function valueByX($x)
	{
		$sum = 0;
		$rev = array_reverse($this->coeffs['a']);
		foreach($rev as $value) {
			$sum*= $x;
			$sum+= $value;
		}
		return $sum;
	}
}