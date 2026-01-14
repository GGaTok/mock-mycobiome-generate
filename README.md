# Mock mycobiome Analysis Pipeline (260113)

분석하고자 하는 genome에서 종 특이적이고 공통서열에 대한 resolution을 확인하기 위해 In silico 에서 인위적인 fastq( Miseq sequencing 결과와 유사한) data를 만드는  pipeline



## 0. 개요

코드는 크게 Genome down ~ merge 까지 5단계로 나눠져 있으며 각각의 Step에서 코드를 어떻게 사용하는 지 확인하며, 원리 ~ 과정에 대해 정리.



## 1. Reference genome down (from ncbi)

1_reference_genome_down.sh을 실행하게 되면 pretreatment_list.txt에 존재하는 scientific name을 읽고 해당 종의(ncbi에 등록된) reference genome을 다운로드를 진행한다.

다음 forward sequence.fasta와 reverse sequence.fasta를 local blast db로 만들고 다운로드한 genome을 query하여 찾고자 하는 서열의 위치를 확인한다.

그 뒤  forward sequence.fasta 가 hit 한 서열에서 부터 +3000bp를 만일 forward가 hit 하지 못하면 reverse sequence의 -3000bp를 extract하여 정리하도록 코드를 작성했다.

```bash
bash 1_reference_genome_down.sh <forward sequence.fasta> <reverse sequence.fasta>
```


|preteatment_list.txt 예시|
|---|
| **Candida albicans**             |
| **Candida parapsilosis**         |
| **Candida tropicalis**           |
| **Candida kefyr**                |
| **Cladosporium cladosporioides** |
| **Saccharomyces cerevisiae**     |

코드가  성공적으로 작동하면 out_fasta dir와 다운 실패한 종에 대한 log가 남아있는 download_failures.txt 가 생성된다.



```bash
cat out_fasta/*.fasta > Mock.fasta
```

다음 step을 위해  위 코드를 실행 뒤 Mock.fasta를 만들어 준다.



## 2. Random samping.py(선택적)

해당 Step은 분석하고자 하는 종이 많을때 random하게 sampling 을 진행하는 코드이다.

```bash
python 2_random_sampling.py Mock.fasta <sampling 하고자 하는 종의 수>
```

위 코드를 돌릴시 Mock_sampled.fasta 파일이 생성된다.



## 3. Trimming

해당 Step은 random sampling된 fasta 파일에서 target하는 서열을 trimming하는 코드이다.

Primer_set.fasta안에 trimming을 진행하고자 하는 서열 2개를 fasta 형식으로 file을 만들면 된다.

추가적으로 `--max-mm num` 인자를 이용하면 mismatch의 개수를 정할 수 있다.

```bash
python 3_trimming_ver2.py Mock_sampled.fasta <Primer_set.fasta>
```

ITS 1 region만 남길시 아래와 같이 파일을 구성해 주면 된다.

ex) Primer_set.fasta

```
>ITS1
TCCGTAGGTGAACCTGCGG
>ITS2
GCTGCGTTCTTCATCGATGC
```

코드가 정삭적으로 작동시 Mock_sampled_trimmed.fasta 파일이 생성된다. 

※ ITS1 , ITS2, ITS 1+2 region을 모두 보고 싶을시 

Mock_sampled_ITS1_trimmed.fasta, Mock_sampled_ITS2_trimmed.fasta, Mock_sampled_ITS12_trimmed.fasta 3개의 파일을 준비한다.



## 4. Generate reads

해당 Step은 fasta파일에 존재하는 서열을 random한 비율로 fastq 형식으로 파일을 만드는 코드이다. 

일반적으로 fasta 파일을 Quality 40 으로 fastq 로 변환시 Qiime2 dada2 denoising step에서 오류가 발생하기 때문에 이를 해결하고자 Miseq sequencing 결과와 비슷한 quality가 나오도록 (길이가 길어질수록 Quality가 떨어지는) 코드를 설계하였다.

```
python 4_generate_reads_v4.py \
  --input_fastas Mock_sampled_ITS1_trimmed.fasta Mock_sampled_ITS2_trimmed.fasta Mock_sampled_ITS12_trimmed.fasta \
  -n 5 \
  -r 10000 \
  -l 300 \
  -c 4 
```

`-n` : 만들 sample 수 

`-r` : fastq 파일의 read 수 

`-l` : 앞에서 또는 뒤에서 몇 bp인지

`-c` : 코어수

코드가 정상적으로 작동하면

Sample_1~n_Composition.txt

Sample_1 ~ n_Mock_sampled_ITS1 ~ 12_trimmed_1 ~ 2.fastq  

다음과 같은 파일들이 생성된다.

## 5. Merge

해당 Step은 위 step에서 만들어진 fastq 파일을 merge 하는 step이다. 

우리가 진행하고자 하는 실험은 너무 특이적이기 때문에 merge 되는 reads 들은 merge를 진행시키고, 이 외의 read는 linked 해야하기 때문에 이를 위해 만들어진 코드이다.

코드의 경우 read를 우선적으로 merge 시키고 merge 되지 않은 부분에서는 강제적을 linked 시키기 때문에 분석 자체에 문제가 생길 수 있다.

```
python 5_merge_ver2.py .
```

작업하는 dir에 존재하는 모든 fastq 파일을 merge한다.

코드가 정상적으로 작동하면 

Sample_1 ~ n_Mock_sampled_ITS1 ~ 12_trimmed_assembled.fastq

가 생성된다.  이렇게 완성된 mock community를 분석에 사용하면 된다.
