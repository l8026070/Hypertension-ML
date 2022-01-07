<pre>
<?php
ini_set('display_errors', 1); ini_set('display_startup_errors', 1); error_reporting(E_ALL & ~E_NOTICE);

function array_key_first($ar){
	foreach($ar as $k =>$v)
		return $k;
}
$N=1;
//$T=array('tr','ts');
$file='405.clean.nodup.omnionly.bin';
$Res=array();
$R=1000;

$Rs=array();
for($i=0;$i<$N;$i++){
	$f=implode("",file("$file.t$i.profile"));
	$af=explode("\n",$f);	
	foreach($af as $j=>$l){
		$al=explode(" ",preg_replace('/\s+/',' ', trim($l)));
		if (isset($al[2]) && $al[2]>0)
			$Rs[$i][]=array($al[2]-1,number_format($al[5],20,".","")*1000000000000);
		//if ($j>500) break;
	}
}
//print_r($Rs);die();
foreach($Rs as $k => $aR){
	$T=array();
	foreach($aR as $ar)
		$T[$ar[0]][]=$ar[1];
	$Min=min($T[0]);
	$Max=max($T[1]);
	$Del=$Max-$Min;
	$Tot=$Diff=array();
	for($i=0;$i<$R;$i++) 
		foreach($T as $p => $Ar) 
			foreach($Ar as $v)
				if ($v>=$i*$Del/$R+$Min)//simple model! 
					$Tot[$p][$i]++;
	//print_r($Tot);die();
	foreach($Tot as $p =>$Ar)
		for($i=0;$i<$R-1;$i++) 
			$Diff[$p][$i]=$Ar[$i]-$Ar[$i+1];
	//print_r($Diff);die();
	$TP=$TN=$FP=$FN=$FPR=$TPR=array();
	for($i=0;$i<$R;$i++){
		for($j=0;$j<$i;$j++){
			$TN[$i]+=$Diff[0][$j];
			$FN[$i]+=$Diff[1][$j];
		}
		for($j=$i;$j<$R;$j++){
			$FP[$i]+=$Diff[0][$j];
			$TP[$i]+=$Diff[1][$j];			
		}
	}
	//print_r($TN);print_r($FN);print_r($FP);print_r($TP);die();
	$RC=array();
	for($i=0;$i<$R;$i++){
		$TPR[$i]=$TP[$i]/($TP[$i]+$FN[$i]+1);
		$FPR[$i]=$FP[$i]/($FP[$i]+$TN[$i]+1);
		$Ri[$i]=sqrt((1-$TPR[$i])*(1-$TPR[$i])+$FPR[$i]*$FPR[$i]);
		$RC[$FPR[$i]*$R]=$TPR[$i]*$R;
	}
	//print_r($TPR);print_r($FPR);print_r($Ri);print_r($RC);die();
	$AUC=0;
	$v=0;
	for($i=0;$i<$R;$i++){
		if (isset($RC[$i])) $v=$RC[$i];
		$AUC+=$v/$R;
	}
	asort($Ri);
	$I=array_key_first($Ri);
	$TR=$I*$Del/$R+$Min;
	//echo $I.":".$Ri[$I]."=>$TR";die();
	$TP=$TN=$FP=$FN=0;
		
	foreach($T as $p => $Ar) foreach($Ar as $v) {
		if ($v>=$TR){
			if ($p) $TP++;
			else $FP++;
		}else{
			if ($p) $FN++;
			else $TN++;
		}
	}
	$Res[$k]=array($TP,$TN,$FP,$FN,$AUC);
	//break;
}
//print_r($Res);die();
$TP=$TN=$FP=$FN=$AUC=0;
foreach($Res as $k =>$ar){
	$TP +=$ar[0]/$N;
	$TN +=$ar[1]/$N;
	$FP +=$ar[2]/$N;
	$FN +=$ar[3]/$N;
	$AUC+=$ar[4]/$N;
}

//echo "$TP,$TN,$FP,$FN";
$RS['sen']=$TP/($TP+$FN);
$RS['spe']=$TN/($TN+$FP);
$RS['pre']=$TP/($TP+$FP);
$RS['acc']=($TP+$TN)/($TP+$TN+$FP+$FN);
$RS['FS1']=2*$TP/(2*$TP+$FP+$FN);
$RS['AUC']=$AUC/$R;
print_r($RS);

?>