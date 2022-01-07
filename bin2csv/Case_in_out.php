<pre>
<?php

ini_set('memory_limit',-1);
set_time_limit(-1);
error_reporting(E_ALL && ~E_NOTICE);
ini_set('display_errors',1);

$file="top10_rs";//Farhood:farhood_rs, GWAS:gwas_rs, Causal:1, farhood+gwas: farhood_gwas_rs, Top10: top10_rs

$afs=explode("\n",implode("",file("$file.list")));
$ar_rs=array();
foreach($afs as $l) $ar_rs[trim($l)]=0;

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
			foreach($ar_rs as $rs=>$i)
				$in[$I][]=intval($fl[$i]);		}
		$I++;
		if ($I>$IN) break;
	}
	closedir($handle);
}
	
$Data=unserialize(implode("",file("../Phenotype/Data.array")));
$i=0;
foreach($aK as $I=> $CD){
	$ar=$Data[$CD];
	if ($ar['Gender']==0){
		$ar['Gender']=1;
		if(mt_rand(0,9)<5) $ar['Gender']++;
	}
	$Gen=2*($ar['Gender']-1);	
	$s=0;
	foreach($ar['Time'] as $J => $Ar){
		if ($Ar['dbp']>=90 || $Ar['sbp']>=140) $s++;
		if ($Ar['ifhtndrg']==1) $s+=count($ar['Time']);
	}
	$H=0;
	if ($s/count($ar['Time'])>=0.5) $H=1;
	$input[$i]=$in[$I];
	$input[$i][]=$Gen;
	$output[$i][0]=$H;
	
	$i++;
}
/*print_r($in);
print_r($aK);
print_r($input);
print_r($output);
die();*/


$astr=array();
foreach($ar_rs as $rs=>$i){
	$astr[0].="$rs,";
	$astr[1].="continuous,";
	$astr[2].=",";
}
$astr[0].="SEX,";
$astr[1].="discrete,";
$astr[2].=",";

for($i=0;$i<count($output[1]);$i++){ 
	$astr[0].="H$i";
	$astr[1].="discrete";
	$astr[2].="meta";
}

foreach($input as $i=>$ar)
	$astr[]=implode(",",$ar).",".implode(",",$output[$i]);

unset($astr[1]);
unset($astr[2]);
$str=implode("\r\n",$astr);
ob_clean();
header('Content-Type:text/plain;charset=UTF-8');
header('Content-Disposition:attachment;filename="data.csv"');
echo $str;

?>