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
		$x=2*mt_rand(0,$max)/$max;
		if ($max>1) $x-=1;
		$v=($op=='rand')?$x:$op;
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
function ArF($ax,$func){
	$az=array();
	foreach($ax as $j=>$x)
		$az[$j]=ActiveF($x,$func);		
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
function MixG($a1,$a2,$e1,$e2){
	$ar=array();
	$E1=$E2=0.5;
	if ($e1+$e2>0 && 0){
		$E1=$e2/($e1+$e2);
		$E2=$e1/($e1+$e2);
	}
	foreach($a1 as $i =>$a){
		if (is_array($a)) foreach($a as $j=>$v1){
			$v2=$a2[$i][$j];
			$ar[$i][$j]=$v1*$E1+$v2*$E2;
		}else{
			$v2=$a2[$i];
			$ar[$i]=$a*$E1+$v2*$E2;
		}
	}
	return $ar;
}
function NormalA($Ar){
	$R=array();
	$Ar=TransP($Ar);
	foreach($Ar as $j=>$ar){
		$mx=max($ar);
		$mn=min($ar);
		foreach($ar as $i =>$v){
			if ($mx==$mn) $R[$i][$j]=$v;
			else $R[$i][$j]=($v-$mn)/($mx-$mn);
		}
	}
	return $R;
}
function CR2($arR){
	global $ar_time;
	$TS=microtime(1);
	$Ar=array();
	foreach($arR as $i =>$PT){
		$t=$PT[0];
		$p=$PT[1];
		foreach($p as $j=>$v){
			$Ar[$j][0][$i]=$t[$j];
			$Ar[$j][1][$i]=$p[$j];
		}
	}
	//print_r($Ar);
	$R2=array();
	foreach($Ar as $j =>$ar){
		$p=$ar[1];
		$t=$ar[0];
		$m=array_sum($t)/count($t);
		$sp=$sm=0;
		foreach($t as $i=>$y){
			$sp+=pow($y-$p[$i],2);
			$sm+=pow($y-$m,2);
		}
		$R2[$j]=1-($sp/$sm);
	}
	$ar_time['CR2']+=microtime(1)-$TS;
	return $R2;
}
function C_AUC($APR){
	global $ar_time;
	$TS=microtime(1);
	$dt=array();
	foreach($APR as $j => $PR){
		if (!isset($PR[0])) $PR[0]=$PR[1];
		//Train
		$au=$kus=array();
		for($k=1;$k<100;$k++){
			$th=$k/100;
			$tp=$fn=$fp=$tn=0;
			foreach($PR[0] as $i=>$v){//train
				list($va,$pd,$F)=$v;
				if ($va>$th && $pd>$th) $tp++;
				if ($va>$th && $pd<=$th) $fn++;
				if ($va<=$th && $pd>$th) $fp++;
				if ($va<=$th && $pd<=$th) $tn++;
			}
			$au[$k]=array($tp,$fn,$fp,$tn,$th);
		}
		$aus=array();
		foreach($au as $k => $av){
			$fpr=$tpr=0;
			list($tp,$fn,$fp,$tn,$th)=$av;
			if ($fp+$tn>0) $fpr=round(100*$fp/($fp+$tn));
			if ($tp+$fn>0) $tpr=round(100*$tp/($tp+$fn));
			$aus[$fpr][$tpr][]=$k;
		}
		//print_r($aus);
		$K=$ux=$mx=0;	
		foreach($aus as $fpr =>$ar) foreach($ar as $tpr =>$an){
			$d=pow($fpr-100,2)+pow($tpr,2);
			if ($d>=$mx){
				$mx=$d;
				list($K)=$aus[$fpr][$tpr];
			}
		}
		for($i=0;$i<100;$i++){
			for($k=0;$k<101;$k++) 
				if (isset($aus[$i][$k]))
					$kus[$i]=$k;
			if(!isset($kus[$i])){
				if ($i==0) $kus[0]=0;
				else $kus[$i]=$kus[$i-1];
			}
		}
		//print_r($kus);
		for($i=0;$i<100;$i++) $ux+=$kus[$i];
		$ux/=10000;
		foreach($kus as $x =>$y)
			$dt[$j][]="{x:".($x/100).",y:".($y/100)."}";
		//Test
		list($tp,$fn,$fp,$tn,$th)=$au[$K];
		$tp=$fn=$fp=$tn=0;
		foreach($PR[1] as $i=>$v){//test
			list($va,$pd)=$v;			
			if ($va>$th && $pd>$th) $tp++;
			if ($va>$th && $pd<=$th) $fn++;
			if ($va<=$th && $pd>$th) $fp++;
			if ($va<=$th && $pd<=$th) $tn++;
		}
		echo "<hr>Report $j:";
		echo "\nTHR = $th";
		echo "\nFPR = ".		($fp/($fp+$tn));
		echo "\nTPR = ".		($tp/($tp+$fn));
		echo "\nPrecision = ".	($tp/($tp+$fp));
		echo "\nRecall = ".		($tp/($tp+$fn));
		echo "\nAccuracy = ".	(($tp+$tn)/($tp+$fp+$tn+$fn));
		echo "\nF-measure = ".	(2*$tp/(2*$tp+$fp+$fn));
		echo "\nAUC = ".		$ux;		
	}
	$ar_time['C_AUC']+=microtime(1)-$TS;
	return $dt;
}

ini_set('memory_limit',-1);
set_time_limit(-1);

$ar_time=array();
$t = -microtime(1);
echo "<pre>";

//// INPUT&OUTPUT /////

//Random
//$N=100;
//$input=FillMat($N,5);
//$output=FillMat($N,2);
//$input=FillMat($N,5,'rand',1);
//$output=FillMat($N,2,'rand',1);

//XOR
//$input=[[0,0],[0,1],[1,0],[1,1]];
//$output=[[0],[1],[1],[0]];

//IRIS
//$f=implode("",file("../NN/data/iris.tab"));
//$af=explode("\n",trim($f));
//foreach($af as $i=>$l){
//	$d=explode("\t",trim($l));
//	$input[$i]=array($d[0],$d[1],$d[2],$d[3]);
//	$output[$i]=array(intval($d[4]=='Iris-setosa'),intval($d[4]=='Iris-versicolor'),intval($d[4]=='Iris-virginica'));
//}

//Vertebral
//$f=implode("",file("../NN/data/vertebral.tab"));
//$af=explode("\n",trim($f));
//foreach($af as $i=>$l){
//	$d=explode("\t",trim($l));
//	$input[$i]=array($d[0],$d[1],$d[2],$d[3],$d[4],$d[5]);
//	//$output[$i]=array(intval($d[6]!='NO'));
//	$output[$i]=array(intval($d[6]=='DH'),intval($d[6]=='SL'));	
//	//$output[$i]=array(intval($d[6]=='DH'),intval($d[6]=='SL'),intval($d[6]=='NO'));
//}

//Z=X+Rm*Y+Ra
//for($i=0;$i<100;$i++){
//	$x=rand(0,10);
//	$y=rand(0,10);
//	$rm=rand(1,2);
//	$ra=rand(0,10)-5;
//	$z=$x+$rm*$y+$ra;
//	$input[$i]=array($x,$y);
//	$output[$i]=array($z);
//}
//////////////////////////////////////////////
if (count($input)<>count($output) || count($output)==0) SE('IN&OUT Error!');

//$input=NormalA($input);$output=NormalA($output);

//Print_M('Input',$input);Print_M('Output',$output);

//Parameter
$actf='x';
$errf='mse';
$H_L=array();
$maxIteration=1000;
$Genration=1000;
$minError=0.00001;
$testS=100/100;
$MuteR=5/100;
$MuteG=5/100;

echo "\nmaxIteration: $maxIteration, Generation: $Genration, testS: $testS, MuteR: $MuteR\n";

///// LAYER /////
$LN[0]=count($input[0]);
foreach($H_L as $nhl)
	$LN[]=$nhl;
$LN[]=count($output[0]);
$L_N=count($LN)-1;
print_r($LN);

///Test && Sample /////
$IN=count($input);
$tn=round($testS*$IN);
$tn=max(2,$tn);
$test=array_rand($input,$tn);
if (!is_array($test)) $test=array($test);

//Initiate
$A_GEN=array();
for($k=0;$k<$Genration;$k++) for($j=0;$j<count($LN)-1;$j++){
	$A_GEN[$k]['W'][$j]=FillMat($LN[$j],$LN[$j+1]);
	$A_GEN[$k]['B'][$j]=FillMat(1,$LN[$j+1]);	
	$A_GEN[$k]['E'][$j]=1;
}
foreach($A_GEN as $k =>$Gen){
	$er=array();
	foreach($input as $i =>$av) if (!in_array($i,$test) || $testS==1){
		list($a,$z)=Forward($av,$Gen['W'],$Gen['B'],$actf);
		$e=Compute($a[$L_N],$output[$i],$errf);
		foreach($e as $j =>$v)
			$er[$j]+=$v/$IN;
	}
	$A_GEN[$k]['E']=$er;
}
//Print_M('A',$A_GEN);

///// RUN ///////
$ers=$es=$Ar_C=array();
$I=0;
$MI=floor($maxIteration/100);
$perr=0;
$MR=$MuteR;
while($I<$maxIteration){
	$I++;
	/// Random ///
	$arN=array();
	for($k=0;$k<$Genration;$k++)
		$arN[$k]=array_rand($A_GEN,2);
	/// Genration ///
	for($k=0;$k<$Genration;$k++){
		for($j=0;$j<count($LN)-1;$j++){
			if ((mt_rand(0,100)/100)<$MR){
				$A_GEN[$k+$Genration]['W'][$j]=FillMat($LN[$j],$LN[$j+1]);
				$A_GEN[$k+$Genration]['B'][$j]=FillMat(1,$LN[$j+1]);
			}else{
				$A_GEN[$k+$Genration]['W'][$j]=MixG($A_GEN[$arN[$k][0]]['W'][$j],$A_GEN[$arN[$k][1]]['W'][$j],$A_GEN[$arN[$k][0]]['E'][$j],$A_GEN[$arN[$k][1]]['E'][$j]);
				$A_GEN[$k+$Genration]['B'][$j]=MixG($A_GEN[$arN[$k][0]]['B'][$j],$A_GEN[$arN[$k][1]]['B'][$j],$A_GEN[$arN[$k][0]]['E'][$j],$A_GEN[$arN[$k][1]]['E'][$j]);;
			}
		}
		/// Error ///
		$Gen=$A_GEN[$k+$Genration];
		$er=array();
		foreach($input as $i =>$av) if (!in_array($i,$test) || $testS==1){
			list($a,$z)=Forward($av,$Gen['W'],$Gen['B'],$actf);
			$e=Compute($a[$L_N],$output[$i],$errf);
			foreach($e as $j =>$v)
				$er[$j]+=$v/$IN;
		}
		$A_GEN[$k+$Genration]['E']=$er;
	}
	foreach($A_GEN as $k =>$Gen){
		$eR[$k]=array_sum($A_GEN[$k]['E']);
	}
	/// Reset ///
	$P_GEN=$A_GEN;
	$A_GEN=array();
	asort($eR);
	//print_r($eR);
	$K=0;
	foreach($eR as $k =>$e) if ($K<$Genration){
		$A_GEN[$K]=$P_GEN[$k];
		$K++;
	}
	
	//$err=array_sum($eR)/count($eR);
	$err=array_sum($A_GEN[0]['E']);
	if ($MI==0 || $I%$MI==0 || $err<$minError || $MR>1) $es[$I]=$err;
	if ($perr==$err) $MR*=(1+$MuteG);
	else $MR=$MuteR;
	$perr=$err;
	$TL=$t+microtime(1);
	echo "$I (".round($TL)."): [".number_format($MR,3)."] $err\n";flush(); @ob_flush();
	if ($err<$minError || $MR>1) break;	
	//print_r($A_GEN);
}
//Print_M('Error',$eR);

/////// MSE //////
//echo "\nI: $I \n";
foreach($es as $x =>$y) $ers[]="{x:$x,y:$y}";
$avgmse=array_sum($es)/count($es);

///// MODEL //////
$W=$A_GEN[0]['W'];
$B=$A_GEN[0]['B'];
//Print_M('A',$A_GEN[0]);

/////// Predict /////
$er=$arR=$PR=array();
foreach($input as $i =>$in){
	$F=intval(in_array($i,$test));
	list($a,$z)=Forward($in,$W,$B,$actf);
	$p=$a[count($z)];
	$arR[$F][$i]=array($output[$i],$p);
	foreach($output[$i] as $j=>$tg){
		$pd=$p[$j];
		$PR[$j][$F][$i]=array($tg,$pd);
	}
	if ($F) $er[$i]=array_sum(Compute($p,$output[$i],$errf));
}
//Print_M('arR',$arR);
//Print_M('Predict',$PR);
//Print_M('Err',$er);
$err=array_sum($er)/count($er);
echo "\nError Predict ($errf): $err\n";

////// R2 & AUC /////
$R2=CR2($arR[1]);
Print_M('R2',$R2);

$dts=array();
$dt=C_AUC($PR);
foreach($dt as $j =>$ad)
	$dts[$j]="{name: \"F_$j\",type: \"line\",showInLegend: true,dataPoints: [".implode(",\n",$ad)."]}";

/////////  TIME /////////
$Tot=$t+microtime(1);
asort($ar_time);
Print_M('Time',$ar_time);
echo "\n Total Time: $Tot";

?>

<!DOCTYPE HTML>
<html>
<head>  
<meta charset="UTF-8">
<script>
window.onload = function () {

var rchart = new CanvasJS.Chart("roc_chart", {
	animationEnabled: true,
	title:{text: "ROC"},
	exportEnabled: true,
	axisY:{maximum: 1,includeZero: true},
	axisX:{maximum: 1,includeZero: true},
	data: [<?php echo implode(",\n",$dts);?>]
});
var echart = new CanvasJS.Chart("mse_chart", {
	animationEnabled: true,
	title:{text: "MSE"},
	exportEnabled: true,
	axisY: {
		gridDashType: "dash",
		stripLines: [{
			value: <?php echo $avgmse;?>,
			label: "Average: <?php echo $avgmse;?>",
			labelFontColor: "#FF0000",
			showOnTop: true,
			labelAlign: "center",
			color: "#FF0000"
		}]
	},

	data: [{type: "line",dataPoints: [<?php echo implode(",\n",$ers);?>]}]
});

rchart.render();
echart.render();
}
</script>
</head>
<body>
<div id="mse_chart" style="height: 300px; max-width: 1000px; margin: 0px auto;"></div>
<div id="roc_chart" style="height: 300px; max-width: 300px; margin: 0px auto;"></div>
<script src="../NN/canvasjs.min.js"></script>
</body>
</html>