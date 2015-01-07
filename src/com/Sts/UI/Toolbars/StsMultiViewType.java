package com.Sts.UI.Toolbars;

import com.Sts.UI.*;

import java.io.*;

/**
 * Created by IntelliJ IDEA.
* User: Tom Lasseter
* Date: Nov 8, 2008
* Time: 11:48:35 PM
* To change this template use File | Settings | File Templates.
*/
public class StsMultiViewType implements Serializable
{
    StsViewItem[] viewItems;
    String name;
    public StsMultiViewType()
    {
    }

    public StsMultiViewType(StsViewItem[] viewItems, String name)
    {
        this.viewItems = viewItems;
        this.name = name;
    }

    public void initializeViewClasses()
    {
        for(StsViewItem viewItem : viewItems)
            viewItem.initializeViewClasses();
    }
}
