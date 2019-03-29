using System;
using System.Reflection;
using System.Collections.Generic;
using System.Linq;
using ilPSP;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.IO
{
    class ReflectionVectorDataSerializer
    {
        readonly IVectorDataSerializer vectorSerializer;

        public ReflectionVectorDataSerializer(IVectorDataSerializer vectorSerializer)
        {
            this.vectorSerializer = vectorSerializer;
        }

        public void SaveVector(object[] vector, Type vectorType, Guid guid)
        {
            MethodInfo method = GetGenericMethod(
                typeof(ReflectionVectorDataSerializer),
                "SaveVector_ValueTypeCapable",
                new[] { vectorType },
                new[] { typeof(IList<object>), typeof(Guid)}
                );
            method.Invoke(this, new object[] { vector, guid});
        }

        public Guid SaveVector(object[] vector, Type vectorType)
        {
            MethodInfo method = GetGenericMethod(
                typeof(ReflectionVectorDataSerializer),
                "SaveVector_ValueTypeCapable",
                new[] { vectorType },
                new[] { typeof(IList<object>)}
                );
            Guid id = (Guid)method.Invoke(this, new object[] { vector});
            return id;
        }

        public Guid SaveVector_ValueTypeCapable<T>(IList<object> vector)
        {
            Guid id = vectorSerializer.SaveVector(vector.Cast<T>().ToArray());
            return id;
        }

        public void SaveVector_ValueTypeCapable<T>(IList<object> vector, Guid id)
        {
            vectorSerializer.SaveVector(vector.Cast<T>().ToArray(), id);
        }

        public IList<object> LoadVector(Guid id, Type vectorType, ref Partitioning part)
        {
            MethodInfo method = GetGenericMethod(
                typeof(ReflectionVectorDataSerializer),
                "LoadVector_ValueTypeCapable",
                new Type[] { vectorType },
                new Type[] { typeof(Guid), typeof( Partitioning) }
                );
            IList<object> data = (IList<object>)method.Invoke(this, new object[] { id, part });
            return data;
        }

        public IList<object> LoadVector_ValueTypeCapable<T>(Guid id, ref Partitioning part)
        {
            IList<T> data = vectorSerializer.LoadVector<T>(id, ref part);
            IList<object> objectData = data.Cast<object>().ToArray();
            return objectData;
        }

        //https://stackoverflow.com/questions/588149/referencing-desired-overloaded-generic-method
        static MethodInfo GetGenericMethod(Type className, string methodName, Type[] genericArgTypes, Type[] argTypes)
        {
            var methods = className.GetMethods();
            MethodInfo method = (from m in className.GetMethods()
                                 where m.Name == methodName &&
                                 m.GetGenericArguments().Length == genericArgTypes.Length &&
                                 ParametersAreEqual(m.GetParameters(), argTypes)
                                 select m
                ).Single();

            MethodInfo genericMethod = method.MakeGenericMethod(genericArgTypes);
            return genericMethod;
        }

        unsafe static bool ParametersAreEqual(ParameterInfo[] genericMethodParameters, Type[] methodParameters)
        {
            if (genericMethodParameters.Length != methodParameters.Length)
            {
                return false;
            }
            for (int i = 0; i < genericMethodParameters.Length; ++i)
            {
                Type genericMethodParameter = genericMethodParameters[i].ParameterType;
                Type methodParameter = methodParameters[i];
                bool equal = true;
                if (!genericMethodParameter.IsGenericType)
                {
                    if (genericMethodParameter.IsByRef)
                    {
                        equal = methodParameter.Equals(genericMethodParameter.GetElementType());
                    }
                    else
                    {
                        equal = methodParameter.Equals(genericMethodParameter);
                    }
                }
                if (!equal)
                {
                    return false;
                }
            }
            return true;
        }
    }
}
