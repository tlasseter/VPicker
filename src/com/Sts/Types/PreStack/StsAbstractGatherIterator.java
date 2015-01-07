package com.Sts.Types.PreStack;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 23, 2008
 * Time: 6:14:17 AM
 * To change this template use File | Settings | File Templates.
 */

/** This iterator begins with the first gather and iterates thru to the last gather.
 *  In 3d this must be inline,crossline ordered sequence.  For 2d, this is a simple
 *  linear sequence from first to last gather in the line.
 */
    abstract public class StsAbstractGatherIterator implements Iterator
    {
        StsDataLineSetFace dataLineSet;
        /** number of samples in each and every non-dead trace */
        int nSamples;

        /** sample data and attribute data for gather (trace-ordered) */
        static StsGatherData gatherData = new StsGatherData();
        /** Iterator interface: returns true if next gather exists */
        abstract public boolean hasNext();
        /** initializes the implemented subclass at start of operation */
        abstract public void initialize();

        StsAbstractGatherIterator(StsDataLineSetFace dataLineSet)
        {
            this.dataLineSet = dataLineSet;
            nSamples = dataLineSet.getNSamples(0);
            initialize();
        }

        /** returns a gatherData structure with data filled in.  Reuse gatherData singleton above in this operation. */
        public StsGatherData next()
        {
            return gatherData;
        }

        public void remove()
        {
        }
    }
