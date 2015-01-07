package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Types.*;
import com.Sts.MVC.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 11, 2010
 * Time: 4:07:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsDynamicPropertyType extends StsPropertyType implements StsSerializable
{
    public StsDynamicPropertyType()
    {
    }

    public StsDynamicPropertyType(StsPropertyTypeDefinition propertyTypeName)
    {
        super(propertyTypeName);
    }

    public StsDynamicPropertyType(StsPropertyTypeDefinition propertyTypeName, boolean persistent)
    {
        super(propertyTypeName, persistent);
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    static public StsDynamicPropertyType constructor(StsPropertyTypeDefinition propertyTypeDefinition)
    {
        return new StsDynamicPropertyType(propertyTypeDefinition);
    }
}
