import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GC123 {
	public static void main(String[] args) {
		String path = "Pathway of your file";
		cleanLineBreak("file name", path);
		readFileByLines("tempFile.fasta", path);
	}

	public static void readFileByLines(String inputfileName, String path) {
		File file = new File(path + inputfileName);
		BufferedReader reader = null;
		BufferedWriter bwgc1 = null;
		BufferedWriter bwgc2 = null;
		BufferedWriter bwgc12 = null;
		BufferedWriter bwgc3 = null;
		try {
			String tempString = null;
			FileWriter writergc1 = new FileWriter(path + "gc1.fasta", true);
			bwgc1 = new BufferedWriter(writergc1);
			FileWriter writergc2 = new FileWriter(path + "gc2.fasta", true);
			bwgc2 = new BufferedWriter(writergc2);
			FileWriter writergc12 = new FileWriter(path + "gc12.fasta", true);
			bwgc12 = new BufferedWriter(writergc12);
			FileWriter writergc3 = new FileWriter(path + "gc3.fasta", true);
			bwgc3 = new BufferedWriter(writergc3);
			reader = new BufferedReader(new FileReader(file));
			while ((tempString = reader.readLine()) != null) {
				float gc1 = 0f;
				int Totail = 0;
				float gc2 = 0f;
				float gc3 = 0f;
				if (tempString.contains(">")) {
					bwgc1.write(tempString + "\t");
					bwgc2.write(tempString + "\t");
					bwgc12.write(tempString + "\t");
					bwgc3.write(tempString + "\t");
				} else {
					int length = tempString.length() / 3;
					for (int j = 0; j < length; j++) {
						String tempCodon = tempString.substring(3 * j,
								3 * j + 3);
						Totail++;
						// gc1中的GC
						if ("C".equals(tempCodon.substring(0, 1))
								|| "G".equals(tempCodon.substring(0, 1))
								|| "c".equals(tempCodon.substring(0, 1))
								|| "g".equals(tempCodon.substring(0, 1))) {
							gc1++;
						}
						// gc2中的GC
						if ("C".equals(tempCodon.substring(1, 2))
								|| "G".equals(tempCodon.substring(1, 2))
								|| "c".equals(tempCodon.substring(1, 2))
								|| "g".equals(tempCodon.substring(1, 2))) {
							gc2++;
						}
						// gc3中的GC
						if ("C".equals(tempCodon.substring(2, 3))
								|| "G".equals(tempCodon.substring(2, 3))
								|| "c".equals(tempCodon.substring(2, 3))
								|| "g".equals(tempCodon.substring(2, 3))) {
							gc3++;
						}
					}
					gc1 = gc1 / Totail;
					gc2 = gc2 / Totail;// gc1,gc2都是不计算(AUG+UGG)，所以
										// gc1Totail=gc2Totail
					gc3 = gc3 / Totail;
					bwgc1.write(String.valueOf(gc1) + "\n");
					bwgc2.write(String.valueOf(gc2) + "\n");
					bwgc12.write(String.valueOf((gc1 + gc2) / 2) + "\n");
					bwgc3.write(String.valueOf(gc3) + "\n");
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				reader.close();
				bwgc3.close();
				bwgc12.close();
				bwgc2.close();
				bwgc1.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public static void cleanLineBreak(String inputFileName, String path) {
		File file = new File(path + inputFileName);
		BufferedReader reader = null;
		FileWriter rwriter = null;
		BufferedWriter bw = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			rwriter = new FileWriter(path + "tempFile.fasta", true);
			bw = new BufferedWriter(rwriter);
			String tempString = null;
			tempString = reader.readLine();
			bw.write(tempString + "\n");
			while ((tempString = reader.readLine()) != null) {
				if (tempString.contains(">")) {
					bw.write("\n" + tempString + "\n");
				} else {
					bw.write(tempString);
				}
			}
			bw.flush();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				bw.close();
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
