package com.Sts.DBTypes;

import com.Sts.Types.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 30, 2007
 * Time: 5:02:20 PM
 * To change this template use File | Settings | File Templates.
 */

/** editedBoundingBox is initialized to original boundingBox and then reduced/extended in size and decimated/subdivided in spacing based on croppedBox.
 *  Specifically, the row, col, and slice range variables in this class are adjusted from the originalBoundingBox to the croppedBoundingBox.
 *  On completion, nRows, nCols, and nCroppedSlices correspond with the editedVolume
 */
public class StsEditedBoundingBox extends StsRotatedGridBoundingBox
{
    public StsRotatedGridBoundingBox originalBoundingBox;

    public StsEditedBoundingBox()
    {
    }

    public StsEditedBoundingBox(boolean persistent)
    {
        super(persistent);
    }

    public StsEditedBoundingBox(StsRotatedGridBoundingBox boundingBox, StsRotatedGridBoundingBox croppedBox, boolean persistent)
    {
        super(persistent);
        this.originalBoundingBox = boundingBox;
        initializeToBoundingBox(croppedBox);
    }
}