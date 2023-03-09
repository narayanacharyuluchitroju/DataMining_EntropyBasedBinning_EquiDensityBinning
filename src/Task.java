public class Task {
    String geneId;
    double r1;
    double r2;
    double split;
    double informationGain;

    public Task(String geneId, double r1, double r2, double split, double informationGain) {
        this.geneId = geneId;
        this.r1 = r1;
        this.r2 = r2;
        this.split = split;
        this.informationGain = informationGain;
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public double getR1() {
        return r1;
    }

    public void setR1(double r1) {
        this.r1 = r1;
    }

    public double getR2() {
        return r2;
    }

    public void setR2(double r2) {
        this.r2 = r2;
    }

    public double getSplit() {
        return split;
    }

    public void setSplit(double split) {
        this.split = split;
    }

    public double getInformationGain() {
        return informationGain;
    }

    public void setInformationGain(double informationGain) {
        this.informationGain = informationGain;
    }
}
