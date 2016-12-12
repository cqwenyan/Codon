package wenyan;

public class StringToArray {
	// 将一条基因变成一个密码子数组
	public String[] genArrays(String genString) {
		int genlength = 0;
		int codonNumber = 0;

		genlength = genString.length();
		codonNumber = genlength / 3;
		String[] genArray = new String[codonNumber];
		for (int i = 0; i < codonNumber; i++) {
			genArray[i] = genString.substring(i * 3, (i + 1) * 3);
		}
		return genArray;
	}

	// 将一条基因变成一个密码子二维数组
	public int[][] twoDiArrays(String[] genArray) {
		int PheCount = 0;
		int LeuCount = 0;
		int SerCount = 0;
		int TyrCount = 0;
		int CysCount = 0;
		int ProCount = 0;
		int HisCount = 0;
		int GinCount = 0;
		int ArgCount = 0;
		int IleCount = 0;
		int ThrCount = 0;
		int AsnCount = 0;
		int LysCount = 0;
		int ValCount = 0;
		int AlaCount = 0;
		int AspCount = 0;
		int GluCount = 0;
		int GlyCount = 0;
		int[][] twoDiArray = new int[19][];
		for (int j = 0; j < genArray.length; j++) {
			if (genArray[j].equals("TTT") || genArray[j].equals("TTC")) {
				PheCount++;
			} else if (genArray[j].equals("TTA") || genArray[j].equals("TTG")
					|| genArray[j].equals("CTT") || genArray[j].equals("CTC")
					|| genArray[j].equals("CTA") || genArray[j].equals("CTG")) {
				LeuCount++;
			} else if (genArray[j].equals("TCT") || genArray[j].equals("TCC")
					|| genArray[j].equals("TCA") || genArray[j].equals("TCG")
					|| genArray[j].equals("AGT") || genArray[j].equals("AGC")) {
				SerCount++;
			} else if (genArray[j].equals("TAT") || genArray[j].equals("TAC")) {
				TyrCount++;
			} else if (genArray[j].equals("TGT") || genArray[j].equals("TGC")) {
				CysCount++;
			} else if (genArray[j].equals("CCT") || genArray[j].equals("CCC")
					|| genArray[j].equals("CCA") || genArray[j].equals("CCG")) {
				ProCount++;
			} else if (genArray[j].equals("CAT") || genArray[j].equals("CAC")) {
				HisCount++;
			} else if (genArray[j].equals("CAA") || genArray[j].equals("CAG")) {
				GinCount++;
			} else if (genArray[j].equals("CGT") || genArray[j].equals("CGC")
					|| genArray[j].equals("CGA") || genArray[j].equals("CGG")
					|| genArray[j].equals("AGA") || genArray[j].equals("AGG")) {
				ArgCount++;
			} else if (genArray[j].equals("ATT") || genArray[j].equals("ATC")
					|| genArray[j].equals("ATA")) {
				IleCount++;
			} else if (genArray[j].equals("ACT") || genArray[j].equals("ACC")
					|| genArray[j].equals("ACA") || genArray[j].equals("ACG")) {
				ThrCount++;
			} else if (genArray[j].equals("AAT") || genArray[j].equals("AAC")) {
				AsnCount++;
			} else if (genArray[j].equals("AAA") || genArray[j].equals("AAG")) {
				LysCount++;
			} else if (genArray[j].equals("GTT") || genArray[j].equals("GTC")
					|| genArray[j].equals("GTA") || genArray[j].equals("GTG")) {
				ValCount++;
			} else if (genArray[j].equals("GCT") || genArray[j].equals("GCC")
					|| genArray[j].equals("GCA") || genArray[j].equals("GCG")) {
				AlaCount++;
			} else if (genArray[j].equals("GAT") || genArray[j].equals("GAC")) {
				AspCount++;
			} else if (genArray[j].equals("GAA") || genArray[j].equals("GAG")) {
				GluCount++;
			} else if (genArray[j].equals("GGT") || genArray[j].equals("GGC")
					|| genArray[j].equals("GGA") || genArray[j].equals("GGG")) {
				GlyCount++;
			}
		}
		// 0给AA顺序留着
		twoDiArray[0] = new int[genArray.length];
		twoDiArray[1] = new int[PheCount];
		twoDiArray[2] = new int[LeuCount];
		twoDiArray[3] = new int[SerCount];
		twoDiArray[4] = new int[TyrCount];
		twoDiArray[5] = new int[CysCount];
		twoDiArray[6] = new int[ProCount];
		twoDiArray[7] = new int[HisCount];
		twoDiArray[8] = new int[GinCount];
		twoDiArray[9] = new int[ArgCount];
		twoDiArray[10] = new int[IleCount];
		twoDiArray[11] = new int[ThrCount];
		twoDiArray[12] = new int[AsnCount];
		twoDiArray[13] = new int[LysCount];
		twoDiArray[14] = new int[ValCount];
		twoDiArray[15] = new int[AlaCount];
		twoDiArray[16] = new int[AspCount];
		twoDiArray[17] = new int[GluCount];
		twoDiArray[18] = new int[GlyCount];

		int PheCountSec = 0;
		int LeuCountSec = 0;
		int SerCountSec = 0;
		int TyrCountSec = 0;
		int CysCountSec = 0;
		int ProCountSec = 0;
		int HisCountSec = 0;
		int GinCountSec = 0;
		int ArgCountSec = 0;
		int IleCountSec = 0;
		int ThrCountSec = 0;
		int AsnCountSec = 0;
		int LysCountSec = 0;
		int ValCountSec = 0;
		int AlaCountSec = 0;
		int AspCountSec = 0;
		int GluCountSec = 0;
		int GlyCountSec = 0;
		for (int j = 0; j < genArray.length; j++) {
			// genArrays[j] = genString.substring(j, j + 3);
			if (genArray[j].equals("TTT") || genArray[j].equals("TTC")) {
				twoDiArray[0][j] = 1;
				twoDiArray[1][PheCountSec++] = j;
			} else if (genArray[j].equals("TTA") || genArray[j].equals("TTG")
					|| genArray[j].equals("CTT") || genArray[j].equals("CTC")
					|| genArray[j].equals("CTA") || genArray[j].equals("CTG")) {
				twoDiArray[0][j] = 2;
				twoDiArray[2][LeuCountSec++] = j;
			} else if (genArray[j].equals("TCT") || genArray[j].equals("TCC")
					|| genArray[j].equals("TCA") || genArray[j].equals("TCG")
					|| genArray[j].equals("AGT") || genArray[j].equals("AGC")) {
				twoDiArray[0][j] = 3;
				twoDiArray[3][SerCountSec++] = j;
			} else if (genArray[j].equals("TAT") || genArray[j].equals("TAC")) {
				twoDiArray[0][j] = 4;
				twoDiArray[4][TyrCountSec++] = j;
			} else if (genArray[j].equals("TGT") || genArray[j].equals("TGC")) {
				twoDiArray[0][j] = 5;
				twoDiArray[5][CysCountSec++] = j;
			} else if (genArray[j].equals("CCT") || genArray[j].equals("CCC")
					|| genArray[j].equals("CCA") || genArray[j].equals("CCG")) {
				twoDiArray[0][j] = 6;
				twoDiArray[6][ProCountSec++] = j;
			} else if (genArray[j].equals("CAT") || genArray[j].equals("CAC")) {
				twoDiArray[0][j] = 7;
				twoDiArray[7][HisCountSec++] = j;
			} else if (genArray[j].equals("CAA") || genArray[j].equals("CAG")) {
				twoDiArray[0][j] = 8;
				twoDiArray[8][GinCountSec++] = j;
			} else if (genArray[j].equals("CGT") || genArray[j].equals("CGC")
					|| genArray[j].equals("CGA") || genArray[j].equals("CGG")
					|| genArray[j].equals("AGA") || genArray[j].equals("AGG")) {
				twoDiArray[0][j] = 9;
				twoDiArray[9][ArgCountSec++] = j;
			} else if (genArray[j].equals("ATT") || genArray[j].equals("ATC")
					|| genArray[j].equals("ATA")) {
				twoDiArray[0][j] = 10;
				twoDiArray[10][IleCountSec++] = j;
			} else if (genArray[j].equals("ACT") || genArray[j].equals("ACC")
					|| genArray[j].equals("ACA") || genArray[j].equals("ACG")) {
				twoDiArray[0][j] = 11;
				twoDiArray[11][ThrCountSec++] = j;
			} else if (genArray[j].equals("AAT") || genArray[j].equals("AAC")) {
				twoDiArray[0][j] = 12;
				twoDiArray[12][AsnCountSec++] = j;
			} else if (genArray[j].equals("AAA") || genArray[j].equals("AAG")) {
				twoDiArray[0][j] = 13;
				twoDiArray[13][LysCountSec++] = j;
			} else if (genArray[j].equals("GTT") || genArray[j].equals("GTC")
					|| genArray[j].equals("GTA") || genArray[j].equals("GTG")) {
				twoDiArray[0][j] = 14;
				twoDiArray[14][ValCountSec++] = j;
			} else if (genArray[j].equals("GCT") || genArray[j].equals("GCC")
					|| genArray[j].equals("GCA") || genArray[j].equals("GCG")) {
				twoDiArray[0][j] = 15;
				twoDiArray[15][AlaCountSec++] = j;
			} else if (genArray[j].equals("GAT") || genArray[j].equals("GAC")) {
				twoDiArray[0][j] = 16;
				twoDiArray[16][AspCountSec++] = j;
			} else if (genArray[j].equals("GAA") || genArray[j].equals("GAG")) {
				twoDiArray[0][j] = 17;
				twoDiArray[17][GluCountSec++] = j;
			} else if (genArray[j].equals("GGT") || genArray[j].equals("GGC")
					|| genArray[j].equals("GGA") || genArray[j].equals("GGG")) {
				twoDiArray[0][j] = 18;
				twoDiArray[18][GlyCountSec++] = j;
			}
		}
		return twoDiArray;
	}
}
