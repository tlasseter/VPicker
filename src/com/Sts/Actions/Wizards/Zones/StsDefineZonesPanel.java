package com.Sts.Actions.Wizards.Zones;

import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.event.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Apr 26, 2010
 * Time: 12:00:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsDefineZonesPanel extends StsJPanel
{
    StsTablePanel zoneTable;

    public StsDefineZonesPanel()
    {
        zoneTable = new StsTablePanel();
        zoneTable.setTitle("Select horizon pairs to define units:");
        zoneTable.addColumns(new Object[] { "Zone", "Top", "Base" });
        zoneTable.setPreferredSize(new Dimension(400,300));
        this.add(zoneTable);
    }
}
