package com.Sts.Types.PreStack;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 21, 2008
 * Time: 12:07:37 PM
 * To change this template use File | Settings | File Templates.
 */

/** This class is used by the gatherIterator implemented by the focus data object to transfer gather data
 *  to the S2S I/O processing class (StsPreStackFocusLine3d or StsPreStackFocusLine2d).
 *  In the gatherIterator it is a singleton (static).
 */
public class StsGatherData
{
    /** number of traces in this gather */
    public int nTraces;
    /** number of attributes in each trace header */
    public int nAttributes;
    /** attributes for traces in this gather; dimensioned [nTraces, nAttributes] */
    public double[] traceOrderedAttributes;
    /** singly dimensioned array of gather traces data in trace-order, i.e., ordered [nTraces, nSamples] */
    public float[] traceData;

    public StsGatherData()
    {
    }

    /** for each new gather, set the number of traces
     *
     * @param nTraces number of traces in this gather
     */
    public void setNTraces(int nTraces)
    {
        this.nTraces = nTraces;
    }

    public int getNTraces() { return nTraces; }

    /** for each new gather, set the trace data
     *
     * @param traceData trace data for this gather in trace order
     */
    public void setTraceData(float[] traceData)
    {
        this.traceData = traceData;
    }

    public float[] getTraceData() { return traceData; }

    /** for each new gather, set the attribute data
     *
     * @param traceOrderedAttributes attribute data for gather in trace order, i.e., [nTraces][nAttributes]
     */
    public void setTraceOrderedAttributes(double[] traceOrderedAttributes)
    {
        this.traceOrderedAttributes = traceOrderedAttributes;
    }

    public double[] getTraceOrderedAttributes() { return traceOrderedAttributes; }

	public int getNAttributes() {
		return nAttributes;
	}

	public void setNAttributes(int attributes) {
		nAttributes = attributes;
	}
}