<pre>
<?php

ini_set('memory_limit',-1);
set_time_limit(-1);
error_reporting(E_ALL && ~E_NOTICE);
ini_set('display_errors',1);


$ar_rs=array();
//$afs=explode("\n",implode("",file("farhood_rs.list"))); //Farhood
//$afs=explode("\n",implode("",file("gwas_rs.list"))); //GWAS
//$afs=explode("\n",implode("",file("1.list"))); //Causal
$afs=explode("\n",implode("",file("farhood_gwas_rs.list"))); //farhood+gwas


foreach($afs as $l) $ar_rs[trim($l)]=0;

$path='../Genome/mm/';
$IN=max(intval($_GET['in']),20);

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
		//foreach($fl as $i =>$v) if ($i>0 && $i<$IF)				
		foreach($ar_rs as $rs=>$i)
			$in[$I][$rs]=intval($fl[$i]);//rec: $v>1;dom: $v>0;add: $v //$i-1
		$I++;
		if ($I>$IN) break;
	}
	closedir($handle);
}

$Data=unserialize(implode("",file("../Phenotype/Data.array")));
$arv=array();
foreach($aK as $I=> $CD){
	$ar=$Data[$CD];
	if ($ar['Gender']==0){
		$ar['Gender']=1;
		if(mt_rand(0,9)<5) $ar['Gender']++;
	}
	$Gen=2*($ar['Gender']-1);
	foreach($ar['Time'] as $J => $Ar){
		if ($Ar['dbp']<90 && $Ar['sbp']>140) continue;//Isolated
		if ($Ar['dbp']>=90 || $Ar['sbp']>=140 || $Ar['ifhtndrg']==1)
			foreach($in[$I] as $rs =>$v)
				$arv[$rs][$J][$v]++;
	}
}

//print_r($arv);die();

$astr=array();
$astr[0]="RS,PHASE,Ref,Het,Hom";
$astr[1]="string,discrete,continuous,continuous,continuous";
$astr[2]="meta,meta,,,";

foreach($ar_rs as $rs=>$i) foreach($arv[$rs] as $J=>$ar){
	for($i=0;$i<3;$i++) $ar[$i]=intval($ar[$i]);
	$astr[]="$rs,$J,{$ar[0]},{$ar[1]},{$ar[2]}";
}

$str=implode("\r\n",$astr);
ob_clean();
header('Content-Type:text/plain;charset=UTF-8');
header('Content-Disposition:attachment;filename="data.csv"');
echo $str;

?>