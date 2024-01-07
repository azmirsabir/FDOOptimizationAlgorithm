public class Main {
    public static void main(String []args){
        FDOSingleObjective fdo = new FDOSingleObjective("FX1");
        fdo.setWeightFactor(0.0);
        fdo.setShouldPrintGlobalBestForEachTurn(true);
        fdo.setAgentRemainInsideBoundery(true);
        fdo.setNumberOfAgents(30);
        fdo.setMaxIteration(500);
        fdo.setTurns(1);
        fdo.runFDO();
        fdo.printMean();
        fdo.printStandardDeviation();
    }
}
