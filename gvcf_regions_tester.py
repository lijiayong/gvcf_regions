c_line1 = '1	10001	.	T	<CGA_CNVWIN>	.	.	NS=1;CGA_WINEND=12000	GT:CGA_GP:CGA_NP:CGA_CP:CGA_PS:CGA_CT:CGA_TS:CGA_CL:CGA_LS	.:0.59:1.26:.:0:.:0:0.995:163'
c_line2 = '1	10001	.	T	<CGA_NOCALL>	.	.	END=11038;NS=1;AN=0	GT:PS	./.:.'
c_line3 = '1	11044	.	AAGTCGCACGGCGCCGGGCTGGGGCGGGGGGAGGGTGGCG	.	.	.	NS=1;AN=0	GT:PS	./.:.'
c_line4 = '1	38232	.	A	G	.	.	NS=1;AN=2;AC=2;CGA_XR=dbsnp.131|rs77823476&dbsnp.86|rs806727;CGA_FI=645520|NR_026818.1|FAM138A|TSS-UPSTREAM|UNKNOWN-INC;CGA_SDO=6	GT:PS:FT:GQ:HQ:EHQ:CGA_CEHQ:GL:CGA_CEGL:DP:AD:CGA_RDP	1/1:.:PASS:64:64,136:64,136:48,45:-136,-64,0:-48,-45,0:3:3,3:0'
f_line1 = 'chr22	14257	.	CTG	C	26.01	PASS	AC=2;AF=1;AN=2;DECOMPOSED;GC=53.47;HRun=0;LEN=2;MQ=58.68;NS=1;QD=8.67;TYPE=del	GT:AO:DP:GQ:PL:QA:QR:RO	1|1:2:2:19:38,6,0:40:0:0'
f_line2 = 'chr22	14273	.	C	<*>	0	PASS	AC=0;AF=nan;AN=0;DP=3;END=14282;GC=52.48;MIN=3;MQ=71.25;NS=1	DP:GQ:MIN:QA:QR	3:99:3:0:1203'
g_line1 = '1	63735	.	CCTA	C,<NON_REF>	41.73	.	BaseQRankSum=0.736;ClippingRankSum=0.736;DP=3;MLEAC=1,0;MLEAF=0.500,0.00;MQ=41.10;MQ0=0;MQRankSum=0.736;ReadPosRankSum=0.736	GT:AD:DP:GQ:PL:SB	0/1:1,2,0:3:45:81,0,45,84,51,135:0,1,0,2'
g_line2 = '1	63739	.	C	<NON_REF>	.	.	END=63739	GT:DP:GQ:MIN_DP:PL	0/0:3:0:3:0,0,0'
p_line1 = '1	30867	.	CCTCT	C	49	PASS	BRF=0.0;FR=0.9990;HP=2;HapScore=2;MGOF=20;MMLQ=16;MQ=60.0;NF=0;NR=1;PP=49;QD=53.0;SC=TCTGTGTCTCCCTCTCTCTCT;SbPval=1.0;Source=Platypus;TC=1;TCF=0;TCR=1;TR=1;WE=30879;WS=30857	GT:GL:GOF:GQ:NR:NV	./.:-5.1,-0.0,0.0:20:3:1:1'
p_line2 = '1	30872	.	C	N	5	REFCALL	BRF=.;END=31871;FR=.;FS=.;HP=.;HapScore=.;MGOF=.;MMLQ=.;MQ=.;NF=.;NR=.;PP=.;QD=.;ReadPosRankSum=.;SC=.;START=.;SbPval=.;Size=1000;Source=.;TC=.;TCF=.;TCR=.;TR=.;WE=.;WS=.	GT:GL:GOF:GQ:NR:NV	./.:-1,-1,-1:-1:-1:2:0'

import unittest
import gvcf_regions

class gvcfRegionsTest(unittest.TestCase):
    
    def test_get_bed_region(self):
        self.assertEqual(gvcf_regions.get_bed_region(c_line2), (10000, 11038))
        self.assertEqual(gvcf_regions.get_bed_region(c_line4), (38231, 38232))
        self.assertEqual(gvcf_regions.get_bed_region(f_line1), (14256, 14259))
        self.assertEqual(gvcf_regions.get_bed_region(g_line2), (63738, 63739))
    
    def test_get_GQ(self):
        self.assertEqual(gvcf_regions.get_GQ(c_line4), 64)
        self.assertEqual(gvcf_regions.get_GQ(f_line1), 19)
        self.assertEqual(gvcf_regions.get_GQ(f_line2), 99)
        self.assertEqual(gvcf_regions.get_GQ(g_line1), 45)
        self.assertEqual(gvcf_regions.get_GQ(g_line2), 0)
        
    def test_is_considered(self):
        self.assertEqual(gvcf_regions.is_considered(c_line1, ['CNV']), False)
        self.assertEqual(gvcf_regions.is_considered(c_line1, ['CNV', 'INS:ME']), False)
        self.assertEqual(gvcf_regions.is_considered(c_line1, None), True)
        self.assertEqual(gvcf_regions.is_considered(c_line2, ['CNV']), True)
        
    def test_is_called(self):
        self.assertEqual(gvcf_regions.is_called(f_line1, 5, None, ['PASS']), True)
        self.assertEqual(gvcf_regions.is_called(g_line1, 45, None, None), True)
        self.assertEqual(gvcf_regions.is_called(g_line1, 46, None, None), False)
        self.assertEqual(gvcf_regions.is_called(g_line1, 45, None, ['PASS']), False)
        self.assertEqual(gvcf_regions.is_called(p_line1, None, 49, ['PASS']), True)
        self.assertEqual(gvcf_regions.is_called(p_line1, None, 50, ['PASS']), False)
        self.assertEqual(gvcf_regions.is_called(p_line1, None, 49, ['REFCALL']), False)
        self.assertEqual(gvcf_regions.is_called(p_line1, None, 49, ['PASS', 'REFCALL']), True)

if __name__ == '__main__':
    unittest.main()
