package com.Sts.Types.PreStack;

import javax.media.opengl.GL;

import com.Sts.DBTypes.StsPreStackLineSetClass;
import com.Sts.DBTypes.StsPreStackMicroseismicSet;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: S2S development</p>
 * <p/>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p/>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

/** A superGather consists of a single gather or a group of gathers in a rectangle or cross configuration.
 *  Row and col are the indices of this superGather.  Computations are done on a gather-by-gather basis,
 *  so gatherRow and gatherCol are the indices of the current gather being computed.  gatherData is the
 *  trace data for the current gather at gatherRow & gatherCol.
 */
public class StsMicroseismicSuperGather extends StsSuperGather
{
	transient StsMicroseismicGatherFile file = null;

    public StsMicroseismicSuperGather()
    {

    }
    
    public StsMicroseismicSuperGather(StsModel model, StsPreStackMicroseismicSet lineSet)
    {
//TJL        super(model, lineSet);
    } 
    
    static public StsMicroseismicSuperGather constructor(StsModel model, StsPreStackMicroseismicSet lineSet, StsMicroseismicGatherFile gFile)
    {
        try
        {
        	StsMicroseismicSuperGather microGather = new StsMicroseismicSuperGather(model, lineSet);
        	microGather.setGatherFile(gFile);
            return microGather;
        }
        catch (Exception e)
        {
            StsException.outputException("StsMicroseismicSuperGather.constructor failed.", e, StsException.WARNING);
            return null;
        }
    } 
    
    public void setGatherFile(StsMicroseismicGatherFile gFile)
    {
    	file = gFile;
    }
    public StsMicroseismicGatherFile getGatherFile() { return file; }
    
    public void initializeGathers()
    {
    	nGathersInSuperGather = 1;
    	gatherType = StsPreStackLineSetClass.SUPER_SINGLE;
        if(gathers == null || gathers.length != nGathersInSuperGather)
        {
            //TJL gathers = new StsMicroseismicGather[0];
            //TJL for(int n = 0; n < nGathersInSuperGather; n++)
            //TJL     gathers[n] = new StsMicroseismicGather(this);
        }
        //TJL initializeGather((StsMicroseismicGather)gathers[0]);
        //computeDrawOffsets();
    }
    
    private void initializeGather(StsMicroseismicGather gather)
    {
        //TJL gather.initialize();
        //TJL gather.setSuperGather(this);
        //TJL centerGather = gather;
        //TJL nSuperGatherTraces = gather.nGatherTraces;
    }   
    
    public void displayWiggleTraces(GL gl, float[][] axisRanges, float pixelsPerXunit, float pixelsPerYunit)
    {
        //TJL centerGather.displayWiggleTraces(gl, axisRanges, pixelsPerXunit, pixelsPerYunit);
    }    
}