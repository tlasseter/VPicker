package com.Sts.DBTypes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsMicroseismicBoxSet extends StsObject
{
    private StsObjectRefList boxes;
    private StsColor stsColor;

    public StsMicroseismicBoxSet()
    {
    }

    public void add(StsBoxSubVolume box)
    {
        if(boxes == null)
        {
            boxes = StsObjectRefList.constructor(4, 4, "boxes", this);
        }
        boxes.add(box);
    }
}
