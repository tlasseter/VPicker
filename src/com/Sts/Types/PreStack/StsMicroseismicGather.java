package com.Sts.Types.PreStack;

import com.Sts.DBTypes.StsPreStackMicroseismicSet;
import com.Sts.UI.StsFilterProperties;
import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsMath;
import com.Sts.Utilities.StsSeismicFilter;
import com.Sts.Utilities.Seismic.*;

/**
 * A superGather consists of a single gather or a group of gathers in a rectangle or cross configuration.
 * Row and col are the indices of this superGather.  Computations are done on a gather-by-gather basis,
 * so gatherRow and gatherCol are the indices of the current gather being computed.  gatherData is the
 * trace data for the current gather at gatherRow & gatherCol.
 */
public class StsMicroseismicGather extends StsGather
{
    public StsMicroseismicGather()
    {
    }

    public StsMicroseismicGather(StsMicroseismicSuperGather superGather)
    {
        super(superGather);
    }
}