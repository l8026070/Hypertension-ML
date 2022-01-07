<?php

include("func.php");

function DropW($i,$j,$k,$v){
	//i: Layer j: input Node  k: output Node
	//if (mt_rand(0,100)<50) return 1;
	return $v;
}
function DropB($i,$j,$v){
	//i: Layer j: Node
	//return 1;
	return $v;
}
function UpdateP($x,$out,$z,$a,$W,$B,$actf,$errf,$LearnRate,$GenRate){
	global $ar_time;
	$TS=microtime(1);
	$Mopt=1;
	//Print_M('x',$x);Print_M('out',$out);Print_M('z',$z);Print_M('a',$a);
	//Print_M('W',$W);Print_M('B',$B);
	$Z=count($z);
	for($i=$Z;$i>0;$i--){
		if ($i==$Z) 
			$da[$i]=Ar2Mat(Compute($a[$Z],$out,$errf,1));
		else 
			$da[$i]=MultM($dz[$i+1],TransP($W[$i]),$Mopt);
		$dz[$i]=MultS($da[$i],Ar2Mat(ArF($z[$i],$actf,1)));
		$dw[$i-1]=MultM(TransP(Ar2Mat($a[$i-1])),$dz[$i],$Mopt);
		$db[$i-1]=$dz[$i];
	}
	//Print_M('da',$da);Print_M('dz',$dz);Print_M('dw',$dw);Print_M('db',$db);
	$GR=$GenRate*$LearnRate;
	foreach($W as $i =>$al) foreach($al as $j=>$aw) foreach($aw as $k=>$w){
		$W[$i][$j][$k]=DropW($i,$j,$k,$w*(1-$GR)-$LearnRate*$dw[$i][$j][$k]);
	}
	foreach($B as $i =>$ab) foreach($ab as $j=>$b){
		$B[$i][$j]=DropB($i,$j,$b*(1-$GR)-$LearnRate*$db[$i][0][$j]);
	}	
	//Print_M('W',$W);Print_M('B',$B);
	$ar_time['UpdateP']+=microtime(1)-$TS;
	return array($W,$B);	
}
function StatA($ar){
	$m=array_sum($ar)/count($ar);
	$s2=0;
	foreach($ar as $j =>$v) 
		$s2+=pow($v-$m,2);
	$s=sqrt($s2/(count($ar)-1));
	return array($m,$s);
}
function NormalA($Ar,$type='minmax'){
	$R=array();
	$Ar=TransP($Ar);
	foreach($Ar as $j=>$ar) switch($type){
		case 'minmax':
			$mx=max($ar);
			$mn=min($ar);
			foreach($ar as $i =>$v)
				$R[$i][$j]=($v-$mn)/($mx-$mn);
		break;
		case 'gaussian':
			list($m,$s)=StatA($ar);
			foreach($ar as $i =>$v)
				$R[$i][$j]=($v-$m)/$s;
		break;
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
		//Train
		if (!isset($PR[0])) $PR[0]=$PR[1];
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

//$x=FillMat(500,400);$y=FillMat(400,600);MultM($x,$y);MultMP($x,$y);print_r($ar_time);die();

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
	//one-hot
//	$output[$i]=array(intval($d[4]=='Iris-setosa'),intval($d[4]=='Iris-versicolor'),intval($d[4]=='Iris-virginica'));
//}

//Vertebral
//$f=implode("",file("../NN/data/vertebral.tab"));
//$af=explode("\n",trim($f));
//foreach($af as $i=>$l){
//	$d=explode("\t",trim($l));
//	$input[$i]=array($d[0],$d[1],$d[2],$d[3],$d[4],$d[5]);
//	//one-hot
//	//$output[$i]=array(intval($d[6]!='NO'));
//	$output[$i]=array(intval($d[6]=='DH'),intval($d[6]=='SL'));
//	//$output[$i]=array(intval($d[6]=='DH'),intval($d[6]=='SL'),intval($d[6]=='NO'));
//}

//Z=X+Rm*Y+Ra
//for($i=0;$i<100;$i++){
//	$x=rand(0,10);
//	$y=rand(0,10);
//	$rm=rand(1,1);
//	$ra=rand(1,2)+5;
//	$z=$x+$rm*$y+$ra;
//	$input[$i]=array($x,$y);
//	$output[$i]=array($z);
//}
//Y=Rm*X+Ra
//for($i=0;$i<100;$i++){
//	$x=rand(0,10);
//	$rm=rand(1,1);
//	$ra=rand(1,1)-5;
//	$y=$rm*$x+$ra;
//	$input[$i]=array($x);
//	$output[$i]=array($y);
//}

//////////////////////////////////////////////
if (count($input)<>count($output) || count($output)==0) SE('IN&OUT Error!');
//$input=NormalA($input);$output=NormalA($output);
//Print_M('Input',$input);Print_M('Output',$output);
//foreach($input as $i =>$ar) echo implode(",",$ar).",".implode(",",$output[$i])."\n";

//Stat IN&OUT
$ins=TransP($input);
$outs=TransP($output);
$Stat=array();
foreach($ins as $i =>$ar)
	$Stat['IN'][$i]=StatA($ar);
foreach($outs as $i =>$ar)
	$Stat['OUT'][$i]=StatA($ar);
//Print_M('Stat',$Stat);

//Parameter
$actf='x';//relu,tanh,sigm,soft,hard,x
$errf='mse';//mse,abs
$minError=0.00001;

$H_L=array();
$maxIteration=10000;
$LearnRate=0.001;
$GenRate=0.0001;
$KL=0.99;

$testS=100/100;
$sampleS=100/100;

echo "\nminError: $minError, maxIteration: $maxIteration
LearnRate: $LearnRate, GenRate: $GenRate
testS: $testS, sampleS: $sampleS\n";

///Test && Sample /////
$IN=count($input);
$tn=round($testS*$IN);
$tn=max(2,$tn);
$test=array_rand($input,$tn);
if (!is_array($test)) $test=array($test);
$sn=round($sampleS*$IN);

echo "\nInput: $IN\nTest: $tn\nSample: $sn\n";

///// LAYER /////
$LN[0]=count($input[0]);
foreach($H_L as $nhl)
	$LN[]=$nhl;
$LN[]=count($output[0]);
$L_N=count($LN)-1;
print_r($LN);

//Initiate
for($j=0;$j<count($LN)-1;$j++){
	$W[$j]=FillMat($LN[$j],$LN[$j+1],0);//Zero Start!
	$B[$j]=FillMat(1,$LN[$j+1],0);//Zero Start!
}
//Print_M('W',$W);Print_M('B',$B);die();
///// RUN ///////
$ers=$es=array();
$I=0;
$errp=INF;
$derrp=INF;
$Tfor=$Tupt=0;
$MI=floor($maxIteration/100);
while($I<$maxIteration){
	$I++;
	$er=array();
	$sample=array_rand($input,$sn);
	foreach($input as $i =>$av) if ((!in_array($i,$test) || $testS==1) && in_array($i,$sample)){
		$T0=microtime(1);
		list($a,$z)=Forward($av,$W,$B,$actf);
		$er[$i]=array_sum(Compute($a[$L_N],$output[$i],$errf));
		$T1=microtime(1);
		$Tfor+=$T1-$T0;
		list($W,$B)=UpdateP($av,$output[$i],$z,$a,$W,$B,$actf,$errf,$LearnRate,$GenRate);
		$T2=microtime(1);
		$Tupt+=$T2-$T1;
		//print_r($a[$L_N]);
	}
	//print_r($er);
	$err=array_sum($er)/count($er);
	if ($MI==0 || $I%$MI==0 || $err<$minError) $es[$I]=$err;
	$TL=$t+microtime(1);
	$derr=$err-$errp;
	$serr="<span style='color:green'>$err</span>";
	if ($derr>=0) $serr="<span style='color:red'>$err</span>";
	$serr.=" | $LearnRate ";
	
	if ($derr>=0 || $derrp/$derr<1/$KL) $LearnRate*=$KL;//
	
	$errp=$err;
	$derrp=$derr;
	if ($I%$MI==0) echo "$I (".round($TL)."): $serr\n";flush(); @ob_flush();
	if ($err<$minError || $LearnRate<$minError/10000) break;	
}

/////// MSE //////
//echo "\nI: $I \n";
foreach($es as $x =>$y) $ers[]="{x:$x,y:$y}";
$avgmse=array_sum($es)/count($es);

///// MODEL //////
//Print_M('W',$W);Print_M('B',$B);

///// IMPORTANCE /////
$IMP=array();
foreach($W[0] as $i=>$ar) foreach($ar as $j => $w)
	$IMP[$j][$i]=$Stat['IN'][$i][1]*abs($w);
foreach($IMP as $j =>$ar) arsort($IMP[$j]); 
Print_M('IMP',$IMP);

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
echo "\n Total Time: $Tot \tTfor: ".round(100*$Tfor/$Tot)."% \t Tupt: ".round(100*$Tupt/$Tot)."%\n";

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