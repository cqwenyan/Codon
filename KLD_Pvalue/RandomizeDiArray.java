package wenyan;

import java.util.Random;

public class RandomizeDiArray {

	// 是否用static，twoDiArray[从1开始氨基酸][第几个该AA]=所在位置，twoDiArray[0][位置] = 氨基酸种类;
	public int[][] theShuffle(int[][] twoDiArray) {
		for (int i = 1; i < 19; i++) {
			int tempLength = twoDiArray[i].length;
			Random random = new Random();
			int counter = 0;
			// System.out.println("twoDiArray[i].length"+twoDiArray[i].length+"  "+" i: "+i);
			if (0 == tempLength) {
				continue;
			}
			while (counter <= 10000) {
				// int temp = twoDiArray[i].length;
				for (int j = 0; j < tempLength; j++) {
					counter++;
					// int temp2 = random.nextInt(twoDiArray[i].length);
					exchange(random.nextInt(tempLength), j, twoDiArray, i);
					// if(100000==counter){
					// break;
					// }
				}
			}
		}
		return twoDiArray;
	}

	private void exchange(int p1, int p2, int[][] twoDiArray, int i) {
		int temp = twoDiArray[i][p1];
		twoDiArray[i][p1] = twoDiArray[i][p2];
		twoDiArray[i][p2] = temp;
	}

	public String theShuffleOutput(int[][] twoDiArray, String[] genArray) {
		// 按照输入序列的codon长度生成一个存放数组，用于存放顺序
		int[] resultInt = new int[genArray.length];

		int PhePin = 0;
		int LeuPin = 0;
		int SerPin = 0;
		int TyrPin = 0;
		int CysPin = 0;
		int ProPin = 0;
		int HisPin = 0;
		int GinPin = 0;
		int ArgPin = 0;
		int IlePin = 0;
		int ThrPin = 0;
		int AsnPin = 0;
		int LysPin = 0;
		int ValPin = 0;
		int AlaPin = 0;
		int AspPin = 0;
		int GluPin = 0;
		int GlyPin = 0;

		for (int i = 0; i < genArray.length; i++) {
			// 原始序列的氨基酸顺序
			switch (twoDiArray[0][i]) {
			case 1:
				// 将洗牌后的第一个该氨基酸位置赋给结果数组
				resultInt[i] = twoDiArray[1][PhePin];
				PhePin++;
				break;
			case 2:
				resultInt[i] = twoDiArray[2][LeuPin];
				LeuPin++;
				break;
			case 3:
				resultInt[i] = twoDiArray[3][SerPin];
				SerPin++;
				break;
			case 4:
				resultInt[i] = twoDiArray[4][TyrPin];
				TyrPin++;
				break;
			case 5:
				resultInt[i] = twoDiArray[5][CysPin];
				CysPin++;
				break;
			case 6:
				resultInt[i] = twoDiArray[6][ProPin];
				ProPin++;
				break;
			case 7:
				resultInt[i] = twoDiArray[7][HisPin];
				HisPin++;
				break;
			case 8:
				resultInt[i] = twoDiArray[8][GinPin];
				GinPin++;
				break;
			case 9:
				resultInt[i] = twoDiArray[9][ArgPin];
				ArgPin++;
				break;
			case 10:
				resultInt[i] = twoDiArray[10][IlePin];
				IlePin++;
				break;
			case 11:
				resultInt[i] = twoDiArray[11][ThrPin];
				ThrPin++;
				break;
			case 12:
				resultInt[i] = twoDiArray[12][AsnPin];
				AsnPin++;
				break;
			case 13:
				resultInt[i] = twoDiArray[13][LysPin];
				LysPin++;
				break;
			case 14:
				resultInt[i] = twoDiArray[14][ValPin];
				ValPin++;
				break;
			case 15:
				resultInt[i] = twoDiArray[15][AlaPin];
				AlaPin++;
				break;
			case 16:
				resultInt[i] = twoDiArray[16][AspPin];
				AspPin++;
				break;
			case 17:
				resultInt[i] = twoDiArray[17][GluPin];
				GluPin++;
				break;
			case 18:
				resultInt[i] = twoDiArray[18][GlyPin];
				GlyPin++;
				break;
			default:
				break;
			}
		}
		StringBuilder result = new StringBuilder();
		for (int i = 0; i < genArray.length; i++) {
			if (resultInt[i] == 0) {
				result.append(genArray[i]);
			} else {
				result.append(genArray[resultInt[i]]);
			}
		}
		return result.toString();
	}
}
