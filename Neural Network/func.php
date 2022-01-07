<?php
function Print_M($tit,$ar){
	echo "<hr> $tit: ";
	print_r($ar);
}
function SE($str){
	echo $str;
	die();
}
function FillMat($r,$c,$op='rand',$max=1000){
	$ar=array();
	for($i=0;$i<$r;$i++) for($j=0;$j<$c;$j++){
		$x=mt_rand(0,$max)/$max;
		if ($max>1) $x-=0.5;
		$v=($op==='rand')?$x:$op;
		if ($r==1) $ar[$j]=$v;
		else $ar[$i][$j]=$v;
	}
	return $ar;
}
function ActiveF($x,$func){
	switch($func){
		case 'sigm': return 1/(1+exp(-$x));
		case 'hard': return ($x<=0)?0:1;
		case 'tanh': return (exp($x)-exp(-$x))/(exp($x)+exp(-$x));
		case 'soft': return $x/(1+abs($x));
		case 'relu': return max(0,$x);
		default: return $x;
	}
}
function DifferF($x,$func){
	switch($func){
		case 'sigm': 
			$s=ActiveF($x,$func);
			return $s*(1-$s);
		case 'hard': return 0;
		case 'tanh': 
			$t=ActiveF($x,$func);
			return 1-pow($t,2);
		case 'soft': return 1/pow(1+abs($x),2);
		case 'relu': return ($x<=0)?0:1;
		default: return 1;
	}
}
function ArF($ax,$func,$op=0){
	$az=array();
	foreach($ax as $j=>$x)
		if ($op==0) $az[$j]=ActiveF($x,$func);
		else $az[$j]=DifferF($x,$func);
	return $az;
}
function Forward($x,$W,$B,$func){
	global $ar_time;
	$TS=microtime(1);
	$a=$z=array();
	$a[0]=$x;
	foreach($W as $i =>$ai){
		foreach($ai as $j =>$aj) foreach($aj as $k =>$w) 
			$z[$i+1][$k]+=$w*$a[$i][$j];
		foreach($B[$i] as $k =>$b) 
			$z[$i+1][$k]+=$b;
		$a[$i+1]=ArF($z[$i+1],$func);
	}
	//print_r($a);
	$ar_time['Forward']+=microtime(1)-$TS;
	return array($a,$z);
}
function Compute($pred,$out,$func,$op=0){
	global $ar_time;
	$TS=microtime(1);
	$er=array();
	foreach($pred as $j => $v){
		$d=$v-$out[$j];
		if ($d==0) $er[$j]=0;
		elseif ($func=='abs'){
			if ($op==0) $er[$j]=abs($d);
			elseif ($d<0) $er[$j]=-1;
			else $er[$j]=1;
		}else{
			if ($op==0) $er[$j]=pow($d,2);
			else $er[$j]=2*$d;
		}
	}
	$ar_time['Compute']+=microtime(1)-$TS;
	return $er;
}
function MultS($a,$b){
	global $ar_time;
	$TS=microtime(1);
	if (count($a[0])<>count($b[0])) SE('MultS Er!');
	$ar=array();
	foreach($a[0] as $I =>$v)
		$ar[0][$I]=$v*$b[0][$I];
	//$ar[0]=array_map(function($x,$y){return $x*$y;},$a[0],$b[0]);
	$ar_time['MultS']+=microtime(1)-$TS;
	return $ar;
}
function MultM($a,$b,$op=1){
	global $ar_time;
	$TS=microtime(1);
	if (count($a[0])<>count($b)) SE('MultM Er!');
	$Ar=array();
	if ($op==1){
		foreach($a as $i =>$ar) foreach($ar as $j =>$v1) foreach($b[$j] as $k =>$v2) 
			$Ar[$i][$k]+=$v1*$v2;
	}else{
		$b_=TransP($b);
		foreach($a as $i =>$ar) foreach($b_ as $j =>$br)
			//$Ar[$i][$j]=array_sum(array_map('bcmul', $ar, $br));
			$Ar[$i][$j]=array_sum(array_map(function($x,$y){return $x*$y;},$ar,$br));
	}
	$ar_time['MultM_'.$op]+=microtime(1)-$TS;
	//echo "\n".count($a)."_".count($a[0])." * ".count($b)."_".count($b[0]).": ".(microtime(1)-$TS);
	return $Ar;
}
function TransP($a){
	global $ar_time;
	$TS=microtime(1);
	$Ar=array();
	foreach($a as $i =>$ar) foreach($ar as $j =>$v)
		$Ar[$j][$i]=$v;
	$ar_time['TransP']+=microtime(1)-$TS;
	return $Ar;
}
function Ar2Mat($ar){
	$m=array();
	foreach($ar as $i =>$v)
		$m[0][$i]=$v;
	return $m;
}
function AddM($a,$b){
	global $ar_time;
	$TS=microtime(1);
	if (count($a[0])<>count($b[0])) SE('AddM Er!');
	$ar=array();
	foreach($a as $i =>$Ar) foreach($Ar as $j =>$v)
		$ar[$i][$j]=$v+$b[$i][$j];
	$ar_time['AddM']+=microtime(1)-$TS;
	return $ar;
}
?>