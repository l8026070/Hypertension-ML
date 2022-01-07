<pre>
<?php

ini_set('memory_limit',-1);
set_time_limit(-1);
error_reporting(E_ALL && ~E_NOTICE);
ini_set('display_errors',1);

function FuzzyAge($age){
	$a=array();
	$a[0]=1-($age-17)/10;
	$a[1]=($age<27)?($age-17)/10:1-($age-27)/10;
	$a[2]=($age<37)?($age-27)/10:1-($age-37)/10;
	$a[3]=($age-37)/10;
	foreach($a as $i=>$v){
		if ($v<0) $a[$i]=0;
		if ($v>1) $a[$i]=1;
	}
	return $a;
}
/*Test Dominant,Recessive,Additive
$astr=array();
for($i=0;$i<100;$i++){
	$s=mt_rand(0,2);
	$d=intval($s>0);
	$r=intval($s>1);
	$a=$s;
	$astr[]="$s\t$d\t$r\t$a";
}
$str=implode("\r\n",$astr);
ob_clean();
header('Content-Type:text/plain;charset=UTF-8');
header('Content-Disposition:attachment;filename="data.tab"');
echo $str;
die();
*/


$file="top10_rs";//Farhood:farhood_rs, GWAS:gwas_rs, Causal:1, farhood+gwas: farhood_gwas_rs, Top10: top10_rs

$afs=explode("\n",implode("",file("$file.list")));
$ar_rs=array();
//foreach($afs as $l) $ar_rs[trim($l)]=0;

$path='../Genome/mm/';
$IN=max(intval($_GET['in']),20);

$armp=array();
$I=0;
$input=$in=$output=array();
$aK=array();
	
if (count($ar_rs)>0){

	$spl = new SplFileObject('../../../../../home/code/Genotype/405.clean.nodup.omnionly.bin.map');

	$i=1;
	while(!$spl->eof()){
		$as=explode("\t",trim($spl->fgets()));
		if (isset($ar_rs[$as[1]])) 
			$ar_rs[$as[1]]=$i;
		$i++;
	}
	foreach($ar_rs as $rs =>$i) if ($i==0) unset($ar_rs[$rs]);
	//print_r($ar_rs);echo count($ar_rs);echo "\n".implode(", ",array_keys($ar_rs));die();
}	
if ($handle = opendir(realpath($path))) {
	while (false !== ($file = readdir($handle))){ 
		$afn=explode(".",$file);
		$ext=$afn[1];
		if ($ext!='tab') continue;
		$aK[$I]=$afn[0];
		$f=$path.$file;
		if (count($ar_rs)>0){
			$fl=explode("\t",trim(implode("",file($f))));
			//foreach($fl as $i =>$v) if ($i>0 && $i<$IF)				
			foreach($ar_rs as $rs=>$i)
				$in[$I][]=intval($fl[$i]);//rec: $v>1;dom: $v>0;add: $v //$i-1
		}
		$I++;
		if ($I>$IN) break;
	}
	closedir($handle);
}
	
$Data=unserialize(implode("",file("../Phenotype/Data.array")));
$i=0;
$arp=array();
$ar_dbp=array(80,90,100,110);
$ar_sbp=array(120,140,160,180);
$PH=array();
foreach($aK as $I=> $CD){
	$ar=$Data[$CD];
	if ($ar['Gender']==0){
		$ar['Gender']=1;
		if(mt_rand(0,9)<5) $ar['Gender']++;
	}
	$Gen=2*($ar['Gender']-1);
	foreach($ar['Time'] as $J => $Ar) if ($J<3 && !isset($PH[$CD]) || 1){
		
		if ($Ar['dbp']<90 && $Ar['sbp']>140) continue;//Isolated
		
		$input[$i]=$in[$I];
		$ag=intval($Ar['age']>18);
		$ag+=intval($Ar['age']>40);
		$input[$i][]=$Ar['age'];//$ag
		//$Ag=FuzzyAge($Ar['age']);
		//foreach($Ag as $v) $input[$i][]=$v; 
		$input[$i][]=$Gen;
		$input[$i][]=$CD;
		$input[$i][]=$J;
		//for($k=0;$k<4;$k++)
		//	$output[$i][$k]=intval($Ar['dbp']>$ar_dbp[$k] || $Ar['sbp']>$ar_sbp[$k]);
		$output[$i][0]=intval($Ar['dbp']>=90 || $Ar['sbp']>=140 || $Ar['ifhtndrg']==1);
		$PH[$CD]=$J;
		//$output[$i][0]=intval($Ar['dbp']>90);
		//$output[$i][1]=intval($Ar['sbp']>140);
		//$output[$i][0]=intval($Ar['sbp']>140 || $Ar['dbp']>90);
		//$PP=$Ar['sbp']-$Ar['dbp'];
		//$MAP=$Ar['dbp']+$PP/3;
		//heart rate: $MAP=$Ar['dbp']+0.01*exp(4.14-40.74/HR)*$PP;
		//MAP:		65<MAP<110
		//Hypo:		sbp<90 or dbp<60
		//normal:	sbp<120+10 and dbp<80+5
		//High:		sbp<140 or dbp<90
		//G1 Hyper:	sbp<160 or dbp<100
		//G2 Hyper:	sbp<180 or dbp<110
		//G3 Hyper:	sbp>180 or dbp>110
		//Isolated:	sbp>140 and dbp<90
		//PP:		PP<0.25*sbp or PP>0.5*sbp		
		
		$i++;		
	}
}
//print_r($aK);
//print_r($input);
//print_r($output);
//$op=trim($_GET['op']);
//if ($op=='') $op='NN';
//include("$op.php");


$astr=array();
foreach($ar_rs as $rs=>$i){
	$astr[0].="$rs,";
	$astr[1].="continuous,";
	$astr[2].=",";
}
//$astr[0].="A0,A1,A2,A3,SEX,PID,PHASE,";
//$astr[1].="continuous,continuous,continuous,continuous,discrete,string,discrete,";
$astr[0].="AGE,SEX,PID,PHASE,";
$astr[1].="continuous,discrete,string,discrete,";
$astr[2].=",,meta,meta,";

for($i=0;$i<count($output[1]);$i++){ 
	$astr[0].="H$i,";
	$astr[1].="discrete,";
	$astr[2].="meta,";
}
$astr[0].="RID";
	$astr[1].="string";
	$astr[2].="meta";

foreach($input as $i=>$ar)
	$astr[]=implode(",",$ar).",".implode(",",$output[$i]).",$i";

$str=implode("\r\n",$astr);
ob_clean();
header('Content-Type:text/plain;charset=UTF-8');
header('Content-Disposition:attachment;filename="data.csv"');
echo $str;

?>