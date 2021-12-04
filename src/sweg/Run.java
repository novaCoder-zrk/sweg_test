package sweg;

import sweg.algorithm.*;

import java.io.IOException;
import java.util.Date;

public class Run {

    private SummaryGraphModule module;

    public static void main(String[] args) throws IOException {
        Date today = new Date();
        System.out.println(today);
        final String inputPath = args[0];
        System.out.println("input_path: " + inputPath);
        final String outputPath = args[1];
        System.out.println("output_path: " + outputPath);
        final String sumMode = args[2];
        System.out.println("summarization_mode: " + sumMode);
        System.out.println();

        final SummaryGraphModule module;

        if (sumMode.compareTo("sweg") == 0) {
            module = new SWeG(false);
        } else {
            System.out.println("Invalid command.");
            return;
        }
        int edgeCount = Common.execute(module, inputPath, "\t");
        Common.writeOutputs(module, "output/" + outputPath, edgeCount);
    }
}
