<pre>
<?php

ini_set('memory_limit',-1);
set_time_limit(-1);
error_reporting(E_ALL && ~E_NOTICE);
ini_set('display_errors',1);

$ar_time=array();
$t = -microtime(1);

$ar_rs=array();
//$afs=explode("\n",implode("",file("farhood_rs.list"))); //Farhood
$afs=explode("\n",implode("",file("gwas_rs.list"))); //GWAS
//$afs=explode("\n",implode("",file("1.list"))); //Causal

foreach($afs as $l) $ar_rs[trim($l)]=0;

$IN=max(intval($_GET['in']),20);
$path='../Genome/mm/';

$spl = new SplFileObject('../../../../../home/code/Genotype/405.clean.nodup.omnionly.bin.map');

$i=1;
while(!$spl->eof()){
	$as=explode("\t",trim($spl->fgets()));
	if (isset($ar_rs[$as[1]])) 
		$ar_rs[$as[1]]=$i;
	$i++;
	//if ($i>$IN) break;
}
//print_r($ar_rs);echo count($ar_rs);echo "\n".implode(", ",array_keys($ar_rs));die();

$armp=array();
$I=0;
$input=$in=$output=array();
$aK=array();
if ($handle = opendir(realpath($path))) {
	while (false !== ($file = readdir($handle))){ 
		$afn=explode(".",$file);
		$ext=$afn[1];
		if ($ext!='tab') continue;
		$aK[$I]=$afn[0];
		$f=$path.$file;
		$fl=explode("\t",trim(implode("",file($f))));
		foreach($ar_rs as $rs=>$i) if ($i>0)
			$in[$I][]=intval($fl[$i]);//rec: $v>1;dom: $v>0;add: $v //$i-1
		$I++;
		if ($I>$IN) break;
	}
	closedir($handle);
}

$Data=unserialize(implode("",file("../Phenotype/Data.array")));
$i=0;
$DBP=90;
$SBP=140;
$D=6;
$arp=array();

foreach($aK as $I=> $CD){
	$ar=$Data[$CD];
	if ($ar['Gender']==0){//impute sex
		$ar['Gender']=1;
		if(mt_rand(0,9)<5) $ar['Gender']++;
	}
	$Gen=$ar['Gender']-1;
	$DR=array();
	foreach($ar['Time'] as $Ar){		
		if ($Ar['dbp']<$DBP && $Ar['sbp']>$SBP) continue;//Isolated
		$l=intval($Ar['dbp']>=$DBP || $Ar['sbp']>=$SBP || $Ar['ifhtndrg']==1);
		$ag=intval($Ar['age']);
		$DR[$ag]=$l;			
	}
	if (count($DR)==$D){
		$j=0;
		foreach($DR as $ag=>$l){
			$input[$j][$i]=$in[$I];
			$input[$j][$i][]=$Gen;
			$input[$j][$i][]=$ag;
			$output[$j][$i]=$l;
			$j++;
		}
		$i++;
	}
}
//print_r($aK);print_r($input);print_r($output);die();

include("func.php");

function RandG($k){
	global $U_,$V_,$W_,$PG,$W;
	$U_[$k]=FillMat($W,$W);
	$V_[$k]=FillMat($W,$W);
	$W_[$k]=FillMat($W,1);
	$PG[$k]=CalG($k);
}
function CalG($G){
	global $U_,$V_,$W_,$input,$output,$W,$N,$D,$func;
	$Y_=$E_=array();
	$H_[0]=FillMat($N,$W,0);
	for($j=0;$j<$D;$j++){
		$H_[$j]=AddM($H_[$j],MultM($input[$j],$U_[$G]));
		$H_[$j+1]=MultM($H_[$j],$V_[$G]);
		$y=MultM($H_[$j],$W_[$G]);
		foreach($y as $i=>$ar){
			$Y_[$j][$i]=ActiveF($ar[0],$func);
			$E_[$j][$i]=$Y_[$j][$i]-$output[$j][$i];
		}
	}
	$E=0;
	foreach($E_ as $j =>$ar)
		foreach($ar as $i =>$e)
			$E+=$e*$e;
	//print_r($E_);echo $E;
	return sqrt($E);
}
function MixG($g,$g1,$g2){
	global $U_,$V_,$W_,$PG,$W;
	$E1=$E2=0.5;
	if ($PG[$g1]+$PG[$g2]>0){
		$E1=1-$PG[$g1]/($PG[$g1]+$PG[$g2]);
		$E2=1-$PG[$g2]/($PG[$g1]+$PG[$g2]);
	}
	for($i=0;$i<$W;$i++) for($j=0;$j<$W;$j++)
		$U_[$g][$i][$j]=$U_[$g1][$i][$j]*$E1+$U_[$g2][$i][$j]*$E2;
	for($i=0;$i<$W;$i++) for($j=0;$j<$W;$j++)
		$V_[$g][$i][$j]=$V_[$g1][$i][$j]*$E1+$V_[$g2][$i][$j]*$E2;
	for($i=0;$i<$W;$i++) 
		$W_[$g][$i][0]=$W_[$g1][$i][0]*$E1+$W_[$g2][$i][0]*$E2;
	
	$PG[$g]=CalG($g);
}


$W=count($input[0][0]);
$N=count($input[0]);


$U_=$V_=$W_=$PG=array();
$G=10;
$GN=10;
$MR=5;
$func='sigm';

for($k=0;$k<$G;$k++)
	RandG($k);	

for($g=0;$g<$GN;$g++){
	//print_r($PG);
	for($k=0;$k<$G;$k++){
		if (mt_rand(0,100)<$MR) RandG($k+$G);
		else MixG($k+$G,mt_rand(0,$G),mt_rand(0,$G));
	}
	asort($PG);
	$PU_=$PV_=$PW_=$PPG=array();
	$i=0;
	foreach($PG as $k=>$p) if ($i<$G){
		$PU_[$i]=$U_[$k];
		$PV_[$i]=$V_[$k];
		$PW_[$i]=$W_[$k];
		$PPG[$i]=$PG[$k];
		$i++;
	}
	$U_=$PU_;
	$V_=$PV_;
	$W_=$PW_;
	$PG=$PPG;
	echo "<br>".$PG[0];
}
print_r($U_[0]);print_r($V_[0]);print_r($W_[0]);


?>