import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class DPP {
    public static Task entropyBasedBinning(String[] arr1,String[] arr2,int gene)
    {
        //Entropy of s
        int s_negative=0;
        int s_positive=0;
        for(int i=0;i<arr1.length;i++)
        {
            if(arr2[i].equals("negative"))
            {
                s_negative+=1;
            }
            else {
                s_positive+=1;
            }
        }
        double entropy_s =
                -((double)s_positive/(arr1.length))*(Math.log(((double)s_positive/(arr1.length))) / Math.log(2))
                        -((double)s_negative/(arr1.length))*(Math.log(((double)s_negative/(arr1.length))) / Math.log(2));
        //System.out.println(entropy_s);
        Arrays.sort(arr1,new Comparator<String>()
        {
            public int compare(String s1, String s2)
            {
                if (Double.parseDouble(s1)<(Double.parseDouble(s2)))
                {
                    return -1;
                }
                else {
                    return 1;
                }
            }
        });
        double r1 = Double.parseDouble(arr1[0]);
        double r2 = Double.parseDouble(arr1[arr1.length-1]);
        double split = (r2-r1)/2;
        LinkedList<String[]> s1 = new LinkedList<>();
        LinkedList<String[]> s2 = new LinkedList<>();
        for(int i=0;i<arr1.length;i++)
        {
            if(Double.parseDouble(arr1[i])<=split)
            {
                s1.add(new String[]{arr1[i],arr2[i]});
            }
            else {
                s2.add(new String[]{arr1[i],arr2[i]});
            }
        }
        //Entropy of s1
        int s1_negative=0;
        int s1_positive=0;
        for(int i=0;i<s1.size();i++)
        {
            if(s1.get(i)[1].equals("negative"))
            {
                s1_negative+=1;
            }
            else {
                s1_positive+=1;
            }
        }
        //System.out.println(s1_negative+" "+s1_positive+" "+(s1.size()));
       double entropy_s1 = -((double)s1_positive/(s1.size()))*(Math.log(((double)s1_positive/(s1.size()))) / Math.log(2))
               -((double)s1_negative/(s1.size()))*(Math.log(((double)s1_negative/(s1.size()))) / Math.log(2));
        //System.out.print("entropy_s1: "+entropy_s1);
        //Entrophy of s2
        int s2_negative=0;
        int s2_positive=0;

        for(int i=0;i<s2.size();i++)
        {
            if(s2.get(i)[1].equals("negative"))
            {
                s2_negative+=1;
            }
            else {
                s2_positive+=1;
            }
        }
        //System.out.println();
        //System.out.print(s2_negative+" "+s2_positive+" "+(s2.size()));
        double entropy_s2= -((double)s2_positive/(s2.size()))*(Math.log(((double)s2_positive/(s2.size()))) / Math.log(2))-
                ((double)s2_negative/(s2.size()))*(Math.log(((double)s2_negative/(s2.size()))) / Math.log(2));
        //System.out.println("entropy_s2: "+entropy_s2);

        //Information gain:
        double informationGain_s1s2 = ((double)(s1.size())/(arr1.length))*entropy_s1+((double)(s2.size())/(arr1.length))*entropy_s2;

        //System.out.println("informationGain_s1s2: "+informationGain_s1s2);

        //Gain(split,s):
        double gain = entropy_s - informationGain_s1s2;
        //System.out.println("gain: "+gain);
        return new Task("g"+gene,r1,r2,split,informationGain_s1s2);
    }
    public static void Binning(int top_k_genes,int m_intervals, String inputFileName) throws IOException {
        Scanner sc = new Scanner(new File(inputFileName+".csv"));
        File task1a = new File("entropyRank.csv");
        File task1b = new File("entropyItemMap.csv");
        File task1c = new File("entropyItemizedData.txt");
        File task2a = new File("equiDensityItemMap.csv");
        File task2b = new File("equiDensityItemizedData.txt");
        FileWriter fileWriter1 = new FileWriter(task1a);
        FileWriter fileWriter2 = new FileWriter(task1b);
        FileWriter fileWriter3 = new FileWriter(task1c);
        FileWriter fileWriter4 = new FileWriter(task2a);
        FileWriter fileWriter5 = new FileWriter(task2b);
        List<String[]> ll = new LinkedList<>();
        sc.useDelimiter("\n");   //sets the delimiter pattern
        while (sc.hasNext())  //returns a boolean value
        {
            ll.add(sc.next().split(","));
        }
        sc.close();  //closes the scanner
        String[] arr1 = new String[ll.size()];
        String[] arr2 = new String[ll.size()];
        HashMap<String, double[]> task1 = new HashMap<>();
        int counter1 = 0;
        for (int j = 0; j < ll.size(); j++) {
            arr2[j] = ll.get(j)[ll.get(j).length - 1];
        }
        for (int i = 0; i < top_k_genes; i++) {
            for (int j = 0; j < ll.size(); j++) {
                arr1[j] = ll.get(j)[i];
            }
            Task a = entropyBasedBinning(arr1, arr2, i);
            task1.put(a.geneId, new double[]{a.split, a.informationGain});
            if(i<top_k_genes) {
                fileWriter2.append(a.geneId + "," + "'-inf'" + "," + a.split + "," + counter1);
                fileWriter2.append("\n");
                counter1++;
                fileWriter2.append(a.geneId + "," + a.split + "," + "'+inf'" + "," + counter1);
                fileWriter2.append("\n");
                counter1++;
            }
        }
        List<Map.Entry<String, double[]>> list = new LinkedList<Map.Entry<String, double[]>>(task1.entrySet());
        Collections.sort(list, new Comparator<Map.Entry<String, double[]>>() {
            @Override
            public int compare(Map.Entry<String, double[]> o1, Map.Entry<String, double[]> o2) {
                if (o1.getValue()[1] < o2.getValue()[1]) {
                    return 0;
                } else {
                    return 1;
                }
            }
        });
        int count=0;
        LinkedList<String> entropyRank = new LinkedList<>();
        List<String[]> sortedGenesActual = new LinkedList<>();
        List<double[]> sortedGenesActual_parsed = new LinkedList<>();
        List<String[]> sortedGenesItemized = new LinkedList<>();
        for (Map.Entry<String, double[]> entry : list) {
            if(count<top_k_genes) {
                fileWriter1.append(entry.getKey() + "," + entry.getValue()[0] + "," + entry.getValue()[1]);
                fileWriter1.append("\n");
            }
            count++;
            String geneKey = entry.getKey();
            String[] sortedGeneArrayActual = new String[ll.size()];
            double[] sortedGeneArrayActual_parsed = new double[ll.size()];
            String[] sortedGeneArrayItemized = new String[ll.size()];
            int i = Integer.parseInt(geneKey.substring(1));
            for (int j = 0; j < ll.size(); j++) {
                sortedGeneArrayActual[j] = ll.get(j)[i];
                sortedGeneArrayActual_parsed[j] = Double.parseDouble(ll.get(j)[i]);
                if (Double.parseDouble(ll.get(j)[i]) <= entry.getValue()[0]) {
                    sortedGeneArrayItemized[j] = "" + i * 2;
                } else {
                    sortedGeneArrayItemized[j] = "" + (i * 2 + 1);
                }
            }
            sortedGenesItemized.add(sortedGeneArrayItemized);
            sortedGenesActual.add(sortedGeneArrayActual);
            sortedGenesActual_parsed.add(sortedGeneArrayActual_parsed);
        }
        String[] classArray = new String[ll.size()];
        for (int j = 0; j < ll.size(); j++) {
            classArray[j] = ll.get(j)[ll.get(j).length - 1];
        }
        sortedGenesItemized.add(classArray);
        sortedGenesActual.add(classArray);
        //System.out.print(sortedGenesItemized.size()+" "+ sortedGenesItemized.get(0).length);
        for (int i = 0; i < sortedGenesItemized.get(0).length; i++) {
            String genesRow = "";
            for (int j = 0; j < top_k_genes; j++) {
                if (j != ll.get(0).length - 1) {
                    genesRow += sortedGenesItemized.get(j)[i] + ",";
                } else {
                    genesRow += sortedGenesItemized.get(j)[i];
                }
            }
            fileWriter3.append(genesRow);
            fileWriter3.append("\n");
        }

        LinkedList<LinkedList<LinkedList<Double>>> itemizedData = new LinkedList<>();
        for (int i = 0; i < sortedGenesActual_parsed.size(); i++) {
            double[] temp = sortedGenesActual_parsed.get(i);
            Arrays.sort(temp);
            int counter2 = m_intervals;
            int init_bin = 0;
            int bin_size_act = temp.length / m_intervals;
            int bin_size = bin_size_act;
            LinkedList<LinkedList<Double>> genes_m_intervals = new LinkedList<>();
            while (counter2 > 0) {
                LinkedList<Double> bins = new LinkedList<>();
                if (counter2 == 1) {
                    bin_size = temp.length;
                }
                for (int k = init_bin; k < bin_size; k++) {
                    bins.add(temp[k]);
                    if (k < temp.length - 1 && k == bin_size - 1 && temp[k] == temp[k + 1]) {
                        bin_size++;
                    }

                }
                init_bin = bin_size;
                bin_size += bin_size_act;
                counter2--;
                genes_m_intervals.add(bins);
            }
            itemizedData.add(genes_m_intervals);
        }
        int interval_val = 0;
        LinkedList<LinkedList<String>> equiDensityItemizedData = new LinkedList<>();
        for(int i=0;i<top_k_genes;i++)
        {
            LinkedList<Double> interval_stack = new LinkedList<>();
            for(int j=0;j<itemizedData.get(i).size();j++)
            {
                double max = itemizedData.get(i).get(j).get(itemizedData.get(i).get(j).size()-1);
                 interval_stack.add(max);
                if(j==0)
                {
                    fileWriter4.append(list.get(i).getKey()+","+"'-inf"+","+max+","+interval_val);

                }
                else if(j==itemizedData.get(i).size()-1)
                {
                    fileWriter4.append(list.get(i).getKey()+","+interval_stack.get(j-1)+","+"'+inf'"+","+interval_val);
                }
                else {
                    fileWriter4.append(list.get(i).getKey()+","+interval_stack.get(j-1)+","+max+","+interval_val);
                }
                interval_val++;
                fileWriter4.append("\n");
            }
            Collections.sort(interval_stack);
            LinkedList<String> equiDensityItemizedDataArray = new LinkedList<>();
            for(int j=0;j<sortedGenesActual.get(i).length;j++)
            {
                String itemizedVal ="";
                boolean flag = false;
                for(int k=0;k<interval_stack.size();k++)
                {
                    if(Double.parseDouble(sortedGenesActual.get(i)[j])<=interval_stack.get(k))
                    {
                        itemizedVal+=((Integer.parseInt(list.get(i).getKey().substring(1)) * m_intervals) + k);
                        flag=true;
                        break;
                    }
                }
                if(flag==false)
                {
                    itemizedVal+=((Integer.parseInt(list.get(i).getKey().substring(1)) * m_intervals) + m_intervals);
                }
                equiDensityItemizedDataArray.add(itemizedVal);
            }

            equiDensityItemizedData.add(equiDensityItemizedDataArray);
        }
        //System.out.print(equiDensityItemizedData.size()+" "+equiDensityItemizedData.get(0).size());
        for (int i = 0; i < equiDensityItemizedData.get(0).size(); i++) {
            String genesRow = "";
            for (int j = 0; j < top_k_genes; j++) {
                if (j != ll.get(0).length - 1) {
                    genesRow += equiDensityItemizedData.get(j).get(i)+ ",";
                } else {
                    genesRow += equiDensityItemizedData.get(j).get(i);
                }
            }
            fileWriter5.append(genesRow);
            fileWriter5.append("\n");
        }
        fileWriter1.close();
        fileWriter2.close();
        fileWriter3.close();
        fileWriter4.close();
        fileWriter5.close();



    }
    public static void main(String[] args) throws IOException {
        if(args.length!=3)
        {
            System.out.println("To Execute use below command on command Line: ");
            System.out.print("java DPP P1InputData 3 4");
        }
        else {
            System.out.println("Execution started...");
            String inputFileName = args[0];
            int top_k_genes = Integer.parseInt(args[1]);
            int m_intervals = Integer.parseInt(args[2]);
            Binning(top_k_genes, m_intervals, inputFileName);
            System.out.println("Task1 & Task2 Execution Completed");
        }
    }
}