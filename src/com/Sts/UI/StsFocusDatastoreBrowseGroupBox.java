package com.Sts.UI;

import com.Sts.UI.DataTransfer.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Sep 1, 2008
 * Time: 3:12:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsFocusDatastoreBrowseGroupBox extends StsDatastoreBrowseGroupBox
{
    public StsFocusDatastoreBrowseGroupBox(Object transferPanel, String transferPanelDatastoreFieldName)
    {
        super("Datastore selector", transferPanel, transferPanelDatastoreFieldName);
    }

    public StsDatastoreFace[] getDatastores() { return null; }
}