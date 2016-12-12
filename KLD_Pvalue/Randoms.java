package wenyan;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

public class Randoms {// 未对氨基酸是否一样、密码子含量是否一样进行检测
	public static void genIterator(String inputFileName, String outputFileName, int lineLength) {
		try {
			FileWriter writer = new FileWriter(outputFileName);
			BufferedWriter bw = new BufferedWriter(writer);
			String geneName = null;
			File file = new File(inputFileName);
			BufferedReader reader = new BufferedReader(new FileReader(file));
			int counter = 0;
			while ((geneName = reader.readLine()) != null) {
				if (!geneName.contains(">")) {
					bw.write(core(geneName));
				} else if (geneName.contains(">") && counter != 0) {
					bw.write("\n" + geneName + "\n");
				} else {
					bw.write(geneName + "\n");
				}
				counter++;
			}
			bw.flush();
			reader.close();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static String core(String inputSequence) {
		StringToArray myStringToArray = new StringToArray();
		RandomizeDiArray myRandomizeDiArray = new RandomizeDiArray();
		// String inputSequence = "CGTTTATTATTACGCAGA";
		// 一条基因生成一个密码子数组
		String[] mygenArrays = myStringToArray.genArrays(inputSequence);
		// 将该密码子数组放入到固定顺序的二维数组
		int[][] mytwoDiArrays = myStringToArray.twoDiArrays(mygenArrays);
		// twoDiArray[1-18][]分别打乱顺序
		int[][] afterShuffle = myRandomizeDiArray.theShuffle(mytwoDiArrays);
		String resultSequence = myRandomizeDiArray.theShuffleOutput(
				afterShuffle, mygenArrays);
		return resultSequence;
	}

	public static int cleanLineBreak(String path, String inputFileName) {
		File file = new File(path + inputFileName);
		int sensor = 0;
		int lineLength = 0;
		try {
			FileWriter rwriter = new FileWriter(path + "tempFile.fasta");
			BufferedWriter bw = new BufferedWriter(rwriter);
			String tempStringUpper = null;
			BufferedReader reader = new BufferedReader(new FileReader(file));

			// find the first line contain ">".
			tempStringUpper = reader.readLine();
			while (tempStringUpper.contains(">") != true) {
				tempStringUpper = reader.readLine();
			}
			bw.write(tempStringUpper + "\n");
			// end

			while ((tempStringUpper = reader.readLine()) != null) {
				if (sensor == 0) {
					lineLength = tempStringUpper.length();
					sensor++;
				}
				if (tempStringUpper.contains(">")) {
					bw.write("\n" + tempStringUpper + "\n");
				} else {
					if (tempStringUpper.length() >= 65534) {
						System.out
								.println("warning,some gene is too long in your FASTA file,and we lost some base of them");
					}
					bw.write(tempStringUpper);
				}
			}
			bw.flush();
			reader.close();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return lineLength;
	}
}
