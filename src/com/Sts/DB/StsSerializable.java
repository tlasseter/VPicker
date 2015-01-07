package com.Sts.DB;

/** A flag interface indicating efficient Sts serialization should be used.
 *  Sts serialization unlike Java serlialization doesn't write out class schema for every instance written.
 *  Instead for each field, the value or array of values are all that are written out.  If the value is a
 *  class, the class is hierarchically written out using Sts serialization. */

public interface StsSerializable
{
}


