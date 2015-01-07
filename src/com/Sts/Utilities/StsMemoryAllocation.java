package com.Sts.Utilities;

import java.text.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */
abstract public class StsMemoryAllocation
{
    public long maxMemory, maxMemoryToUse;

    double memAllocFraction = 0.8; // fraction of maxMemory to use for seismic cube processing
    double memFreeFraction = 0.00; // fraction of maxMemory which should be free during processing
    boolean debug = false;
    public long nTotalBytes = 0; // amount of memory actually used

    static final private int maxAllowableMemoryMB = 100; // max memory to be used in volume processing operations
    static final private int maxMapBufferSizeMB = 50;    // max size of any mappedByteBuffer in processing or display operations

    static final public long maxAllowableMemory = maxAllowableMemoryMB*1000000;
    static final public long maxMapBufferSize = maxMapBufferSizeMB*1000000;

    abstract void freeMemory();

    protected void initializeMemoryAllocation()
    {
        maxMemory = Runtime.getRuntime().maxMemory();
        maxMemoryToUse = (long)(memAllocFraction*maxMemory);
        maxMemoryToUse = Math.min(maxMemoryToUse, maxAllowableMemory); // jbw this is essential
        System.out.println("Max memory to use: " + maxMemoryToUse);
        checkMemoryStatus("Initial memory status:");
    }

    protected boolean freeMemoryOK()
    {
        long freeMemory = Runtime.getRuntime().freeMemory();
        return freeMemory > (long)(memFreeFraction*maxMemory);
    }

    void checkMemoryStatus(String message)
    {
        NumberFormat numberFormat = NumberFormat.getInstance();
        String freeMemory = numberFormat.format(Runtime.getRuntime().freeMemory());
        String totalMemory = numberFormat.format(Runtime.getRuntime().totalMemory());
        String maxMemory = numberFormat.format(Runtime.getRuntime().maxMemory());

        System.out.println(message + ". freeMemory: " + freeMemory +
            " totalMemory: " + totalMemory + " maxMemory: " + maxMemory);
    }

    public void setMaxMemoryToUse(int mbMemory)
    {
        maxMemoryToUse = (long)mbMemory*1000000;
    }
    public int getMaxMemoryToUse() { return (int)(maxMemoryToUse/1000000); }
}
