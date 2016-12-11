package wenyan;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class KLD {
	public static String[] code_table = { "TTT", "TCT", "TTC", "TCC", "TTA",
			"TCA", "TTG", "TCG", "CTT", "CCT", "CTC", "CCC", "CTA", "CCA",
			"CTG", "CCG", "ATT", "ACT", "ATC", "ACC", "ATA", "ACA", "ATG",
			"ACG", "GTT", "GCT", "GTC", "GCC", "GTA", "GCA", "GTG", "GCG",
			"TAT", "TGT", "TAC", "TGC", "TAA", "TGA", "TAG", "TGG", "CAT",
			"CGT", "CAC", "CGC", "CAA", "CGA", "CAG", "CGG", "AAT", "AGT",
			"AAC", "AGC", "AAA", "AGA", "AAG", "AGG", "GAT", "GGT", "GAC",
			"GGC", "GAA", "GGA", "GAG", "GGG" };

	public static Map<String, String> synAminoAcidBox() {
		HashMap<String, String> synAminoAcidTable = new HashMap<String, String>();
		//除去起始ATG，除去终止TAA,TAG,TGA，除去单个TGG,Trp
		synAminoAcidTable.put("TTT", "Phe");
		synAminoAcidTable.put("TTC", "Phe");
		synAminoAcidTable.put("TCT", "Ser");
		synAminoAcidTable.put("TCC", "Ser");
		synAminoAcidTable.put("TCA", "Ser");
		synAminoAcidTable.put("TCG", "Ser");
		synAminoAcidTable.put("TAT", "Tyr");
		synAminoAcidTable.put("TAC", "Tyr");
		synAminoAcidTable.put("TGT", "Cys");
		synAminoAcidTable.put("TGC", "Cys");
		synAminoAcidTable.put("TTA", "Leu");
		synAminoAcidTable.put("TTG", "Leu");
		synAminoAcidTable.put("CTT", "Leu");
		synAminoAcidTable.put("CTC", "Leu");
		synAminoAcidTable.put("CTA", "Leu");
		synAminoAcidTable.put("CTG", "Leu");
		synAminoAcidTable.put("CCT", "Pro");
		synAminoAcidTable.put("CCC", "Pro");
		synAminoAcidTable.put("CCA", "Pro");
		synAminoAcidTable.put("CCG", "Pro");
		synAminoAcidTable.put("CAT", "His");
		synAminoAcidTable.put("CAC", "His");
		synAminoAcidTable.put("CAA", "Gin");
		synAminoAcidTable.put("CAG", "Gin");
		synAminoAcidTable.put("CGT", "Arg");
		synAminoAcidTable.put("CGC", "Arg");
		synAminoAcidTable.put("CGA", "Arg");
		synAminoAcidTable.put("CGG", "Arg");
		synAminoAcidTable.put("ATT", "Ile");
		synAminoAcidTable.put("ATC", "Ile");
		synAminoAcidTable.put("ATA", "Ile");
		synAminoAcidTable.put("ACT", "Thr");
		synAminoAcidTable.put("ACC", "Thr");
		synAminoAcidTable.put("ACA", "Thr");
		synAminoAcidTable.put("ACG", "Thr");
		synAminoAcidTable.put("AAT", "Asn");
		synAminoAcidTable.put("AAC", "Asn");
		synAminoAcidTable.put("AAA", "Lys");
		synAminoAcidTable.put("AAG", "Lys");
		synAminoAcidTable.put("AGT", "Ser");
		synAminoAcidTable.put("AGC", "Ser");
		synAminoAcidTable.put("AGA", "Arg");
		synAminoAcidTable.put("AGG", "Arg");
		synAminoAcidTable.put("GTT", "Val");
		synAminoAcidTable.put("GTC", "Val");
		synAminoAcidTable.put("GTA", "Val");
		synAminoAcidTable.put("GTG", "Val");
		synAminoAcidTable.put("GCT", "Ala");
		synAminoAcidTable.put("GCC", "Ala");
		synAminoAcidTable.put("GCA", "Ala");
		synAminoAcidTable.put("GCG", "Ala");
		synAminoAcidTable.put("GAT", "Asp");
		synAminoAcidTable.put("GAC", "Asp");
		synAminoAcidTable.put("GAA", "Glu");
		synAminoAcidTable.put("GAG", "Glu");
		synAminoAcidTable.put("GGT", "Gly");
		synAminoAcidTable.put("GGC", "Gly");
		synAminoAcidTable.put("GGA", "Gly");
		synAminoAcidTable.put("GGG", "Gly");
		return synAminoAcidTable;
	}

	// 待用
	/**
	 * public static Map<String,String> allAminoAcidBox(){ HashMap<String,
	 * String> allAminoAcidTable = new HashMap<String, String>();
	 * allAminoAcidTable.put("TGG","Trp"); allAminoAcidTable.put("AUG","Met");
	 * 
	 * allAminoAcidTable.put("TTT","Phe"); allAminoAcidTable.put("TTC","Phe");
	 * allAminoAcidTable.put("TCT","Ser"); allAminoAcidTable.put("TCC","Ser");
	 * allAminoAcidTable.put("TCA","Ser"); allAminoAcidTable.put("TCG","Ser");
	 * allAminoAcidTable.put("TAT","Tyr"); allAminoAcidTable.put("TAC","Tyr");
	 * allAminoAcidTable.put("TGT","Cys"); allAminoAcidTable.put("TGC","Cys");
	 * allAminoAcidTable.put("TTA","Leu"); allAminoAcidTable.put("TTG","Leu");
	 * allAminoAcidTable.put("CTT","Leu"); allAminoAcidTable.put("CTC","Leu");
	 * allAminoAcidTable.put("CTA","Leu"); allAminoAcidTable.put("CTG","Leu");
	 * allAminoAcidTable.put("CCT","Pro"); allAminoAcidTable.put("CCC","Pro");
	 * allAminoAcidTable.put("CCA","Pro"); allAminoAcidTable.put("CCG","Pro");
	 * allAminoAcidTable.put("CAT","His"); allAminoAcidTable.put("CAC","His");
	 * allAminoAcidTable.put("CAA","Gin"); allAminoAcidTable.put("CAG","Gin");
	 * allAminoAcidTable.put("CGT","Arg"); allAminoAcidTable.put("CGC","Arg");
	 * allAminoAcidTable.put("CGA","Arg"); allAminoAcidTable.put("CGG","Arg");
	 * allAminoAcidTable.put("ATT","Ile"); allAminoAcidTable.put("ATC","Ile");
	 * allAminoAcidTable.put("ATA","Ile"); allAminoAcidTable.put("ACT","Thr");
	 * allAminoAcidTable.put("ACC","Thr"); allAminoAcidTable.put("ACA","Thr");
	 * allAminoAcidTable.put("ACG","Thr"); allAminoAcidTable.put("AAT","Asn");
	 * allAminoAcidTable.put("AAC","Asn"); allAminoAcidTable.put("AAA","Lys");
	 * allAminoAcidTable.put("AAG","Lys"); allAminoAcidTable.put("AGT","Ser");
	 * allAminoAcidTable.put("AGC","Ser"); allAminoAcidTable.put("AGA","Arg");
	 * allAminoAcidTable.put("AGG","Arg"); allAminoAcidTable.put("GTT","Val");
	 * allAminoAcidTable.put("GTC","Val"); allAminoAcidTable.put("GTA","Val");
	 * allAminoAcidTable.put("GTG","Val"); allAminoAcidTable.put("GCT","Ala");
	 * allAminoAcidTable.put("GCC","Ala"); allAminoAcidTable.put("GCA","Ala");
	 * allAminoAcidTable.put("GCG","Ala"); allAminoAcidTable.put("GAT","Asp");
	 * allAminoAcidTable.put("GAC","Asp"); allAminoAcidTable.put("GAA","Glu");
	 * allAminoAcidTable.put("GAG","Glu"); allAminoAcidTable.put("GGT","Gly");
	 * allAminoAcidTable.put("GGC","Gly"); allAminoAcidTable.put("GGA","Gly");
	 * allAminoAcidTable.put("GGG","Gly"); return allAminoAcidTable; }
	 **/

	public static void main(String[] args) {
		// Scanner sc = new Scanner(System.in);
		// System.out.println("Please input file name ,e.g.mulberry.fasta");
		// String inputFileName = sc.nextLine().trim();
		// sc.close();
		String path = "C:\\Users\\Administrator\\Desktop\\testRandom\\kld\\";
		noSpaceTempFile("mulberry_300.fasta",path);
		int i = 1;

		// String[] code_table = { "TTT", "ATG", "TAA", "TAG","TGA","GCC"};

		
		for (int q = 0; q < code_table.length; q++) {
			if ((!code_table[q].equals("ATG"))
					&& (!code_table[q].equals("TAA"))
					&& (!code_table[q].equals("TGA"))
					&& (!code_table[q].equals("TAG"))) {
				String iterationCodon = code_table[q];
				double partA = coreKLDPartA(iterationCodon, synAminoAcidBox(),path);
				double partB = coreKLDPartB(iterationCodon, synAminoAcidBox(),path);
				double partAdivB = partA / partB;
				//原文中的概率被均一化为1000个bases碱基，出现的次数
				double partAB = Math.pow(partA / partB, 333.3*partA);
				double x = Math.log(partAB);

				System.out.println("++++" + iterationCodon + "++++");
				System.out.println("partA:" + partA);
				System.out.println("partB:" + partB);
				System.out.println("partAdivB:" + partAdivB);
				System.out.println("partAB:" + partAB);
				System.out.println("log(partAB)" + x);

			}
		}

		// double partA = coreKLDPartA("ATG",synAminoAcidBox());
		// double partB = coreKLDPartB("ATG",synAminoAcidBox());
		// double partAdivB = partA/partB;
		// double partAB = Math.pow(partA / partB, partA);
		// double x = Math.log(partAB);
		//
		// System.out.println("++++"+"ACA"+"++++");
		// System.out.println("partA:"+partA);
		// System.out.println("partB:"+partB);
		// System.out.println("partAdivB:"+partAdivB);
		// System.out.println("partAB:"+partAB);
		// System.out.println("log(partAB)"+x);

	}

	// 计算fasta文件中的基因条数
	private static int countGeneInFasta(String inputFileName) {
		int genNumberInFasta = 0;
		File file = new File(inputFileName);
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString;
			while ((tempString = reader.readLine()) != null) {
				if (tempString.contains(">")) {
					genNumberInFasta++;
				}
			}
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return genNumberInFasta;
	}

	// 生成无换行临时文件tempNoSpace.fasta
	private static void noSpaceTempFile(String inputFileName,String path) {

		File file = new File(path+inputFileName);
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			FileWriter rwriter = new FileWriter(
					path+"tempNoSpace.fasta");
			BufferedWriter bw = new BufferedWriter(rwriter);
			String tempString;
			tempString = reader.readLine();
			bw.write(tempString + "\n");
			while ((tempString = reader.readLine()) != null) {
				if (tempString.contains(">")) {
					bw.write("\n" + tempString + "\n");
				} else {
					bw.write(tempString);
				}
			}
			bw.close();
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	// 将基因组中对应位置的密码子生成string数组，position代表第几个密码子,whichGene为从上到下各个基因。与tempNoSpace.fasta偶联
	private static String[] getPositionCodon(int position,String path) {
		int genNumber = countGeneInFasta(path+"tempNoSpace.fasta");
		File file = new File(
				path+"tempNoSpace.fasta");
		String[] myPositionCodonArray = new String[genNumber];
		String tempString;
		BufferedReader reader =null;
		try {
			reader = new BufferedReader(new FileReader(file));
			int whichGene = 0;
			while ((tempString = reader.readLine()) != null) {
				if (!tempString.contains(">")) {
					myPositionCodonArray[whichGene] = tempString.substring(
							position, position + 3);
					whichGene++;
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}finally{
			try {
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		return myPositionCodonArray;
	}

	// KLD计算的分子。写死为每个bin10个密码，分别5个bin
	private static double coreKLDPartA(String codon,
			Map<String, String> synAminoAcidTable,String path) {

		int theCodonCounter = 0;
		int synCodonCounter = 0;
		double partA = 0;

		// int totalCodonCounter = 0;

		// synCodonCounter统计出与传入codon同义的codon个数（k）
		for (Map.Entry<String, String> entry : synAminoAcidTable.entrySet()) {
			// 比较氨基酸
			if (synAminoAcidTable.get(codon).equals(entry.getValue())) {
				String synCodon = entry.getKey();
				for (int i = 1; i <= 10; i++) {
					String[] temp = getPositionCodon(i,path);
					for (int j = 0; j < temp.length; j++) {
						// totalCodonCounter++;
						if (temp[j].equals(synCodon)) {
							synCodonCounter++;
						}
					}
				}
			}
		}

		// 统计出该种codon的个数（k）
		for (int i = 1; i <= 10; i++) {
			String[] temp = getPositionCodon(i,path);
			for (int j = 0; j < temp.length; j++) {
				// totalCodonCounter++;
				if (temp[j].equals(codon)) {
					theCodonCounter++;
				}
			}
		}
		System.out
				.println("+++" + codon + "synCodonCounter:" + synCodonCounter);
		System.out
				.println("+++" + codon + "theCodonCounter:" + theCodonCounter);
		partA = (((double) theCodonCounter) / synCodonCounter);
		return partA;
	}

	// KLD计算的分母。与tempNoSpace.fasta强偶联
	private static double coreKLDPartB(String codon,
			Map<String, String> synAminoAcidTable,String path) {
		int oneGenomeCodon = 0;
		int synOneGenomeCodon = 0;
		double partB = 0;
		
		BufferedReader reader = null;
		File file = new File(
				path+"tempNoSpace.fasta");
		try {
			String tempString;
			reader = new BufferedReader(new FileReader(file));
			while ((tempString = reader.readLine()) != null) {
				for (int i = 0; i < tempString.length() / 3; i = i + 3) {
					// 该密码子数量
					if (tempString.substring(i, i + 3).equals(codon)) {
						oneGenomeCodon++;
					}
					// 同义密码子数目
					for (Map.Entry<String, String> entry : synAminoAcidTable
							.entrySet()) {
						// 比较氨基酸
						if (synAminoAcidTable.get(codon).equals(
								entry.getValue())) {
							String synCodon = entry.getKey();
							if (tempString.substring(i, i + 3).equals(synCodon)) {
								synOneGenomeCodon++;
							}
						}
					}
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}finally{
			try {
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println("+++" + codon + "oneGenomeCodon:" + oneGenomeCodon);
		System.out.println("+++" + codon + "synOneGenomeCodon:"
				+ synOneGenomeCodon);
		partB = (((double) oneGenomeCodon) / synOneGenomeCodon);

		return partB;
	}

}
