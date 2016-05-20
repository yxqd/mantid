#include "MantidGeometry/Crystal/SpaceGroupRegistration.h"

namespace Mantid {
namespace Geometry {
namespace SpaceGroupRegistration {

SpaceGroupFactoryImpl &factory = SpaceGroupFactory::Instance();

GeneratedSpaceGroupSubscriber GeneratedSpaceGroup(factory);
TransformedSpaceGroupSubscriber TransformedSpaceGroup(factory);
TabulatedSpaceGroupSubscriber TabulatedSpaceGroup(factory);
OrthorhombicSpaceGroupSubscriber OrthorhombicSpaceGroup(factory);
TransformedOrthorhombicSpaceGroupSubscriber
    TransformedOrthorhombicSpaceGroup(factory);

/* Space groups according to International Tables for Crystallography,
 * using the generators specified there.
 */
Subscribe triclinic((TabulatedSpaceGroup(1, "P 1", "x,y,z"),
                     GeneratedSpaceGroup(2, "P -1", "-x,-y,-z")));

/* Monoclinic space groups.
 *
 * Unique axes b and c are given, as well as 3 cell choices where
 * applicable. Since there are only so few monoclinic space group
 * types, all the transformations are given explicitly. It would
 * actually be shorter to provide all group definition explicitly
 * with their generators, but this way the relationship between
 * the groups is more pronounced.
 */
Subscribe monoclinic((
    GeneratedSpaceGroup(3, "P 1 2 1", "-x,y,-z", "P 2"),
    TransformedSpaceGroup(3, "P 1 1 2", "P 1 2 1 | y,z,x"),

    GeneratedSpaceGroup(4, "P 1 21 1", "-x,y+1/2,-z", "P 21"),
    TransformedSpaceGroup(4, "P 1 1 21", "P 1 21 1 | y,z,x"),

    GeneratedSpaceGroup(5, "C 1 2 1", "-x,y,-z", "C 2"),
    TransformedSpaceGroup(5, "A 1 2 1", "C 1 2 1 | -x+z,y,-x", "A 2"),
    TransformedSpaceGroup(5, "I 1 2 1", "C 1 2 1 | -z,y,x-z", "I 2"),
    TransformedSpaceGroup(5, "A 1 1 2", "C 1 2 1 | y,z,x"),
    TransformedSpaceGroup(5, "B 1 1 2", "A 1 1 2 | -y,x-y,z"),
    TransformedSpaceGroup(5, "I 1 1 2", "A 1 1 2 | -x+y,-x,z"),

    GeneratedSpaceGroup(6, "P 1 m 1", "x,-y,z", "P m"),
    TransformedSpaceGroup(6, "P 1 1 m", "P 1 m 1 | y,z,x"),

    GeneratedSpaceGroup(7, "P 1 c 1", "x,-y,z+1/2", "P c"),
    TransformedSpaceGroup(7, "P 1 n 1", "P 1 c 1 | -x+z,y,-x", "P n"),
    TransformedSpaceGroup(7, "P 1 a 1", "P 1 c 1 | -z,y,x-z", "P a"),
    TransformedSpaceGroup(7, "P 1 1 a", "P 1 c 1 | y,z,x"),
    TransformedSpaceGroup(7, "P 1 1 n", "P 1 1 a | -y,x-y,z"),
    TransformedSpaceGroup(7, "P 1 1 b", "P 1 1 a | -x+y,-x,z"),

    GeneratedSpaceGroup(8, "C 1 m 1", "x,-y,z", "C m"),
    TransformedSpaceGroup(8, "A 1 m 1", "C 1 m 1 | -x+z,y,-x", "A m"),
    TransformedSpaceGroup(8, "I 1 m 1", "C 1 m 1 | -z,y,x-z", "I m"),
    TransformedSpaceGroup(8, "A 1 1 m", "C 1 m 1 | y,z,x"),
    TransformedSpaceGroup(8, "B 1 1 m", "A 1 1 m | -y,x-y,z"),
    TransformedSpaceGroup(8, "I 1 1 m", "A 1 1 m | -x+y,-x,z"),

    GeneratedSpaceGroup(9, "C 1 c 1", "x,-y,z+1/2", "C c"),
    TransformedSpaceGroup(9, "A 1 n 1", "C 1 c 1 | -x+z,y,-x", "A n"),
    TransformedSpaceGroup(9, "I 1 a 1", "C 1 c 1 | -z,y,x-z", "I a"),
    TransformedSpaceGroup(9, "A 1 1 a", "C 1 c 1 | y,z,x"),
    TransformedSpaceGroup(9, "B 1 1 n", "A 1 1 a | -y,x-y,z"),
    TransformedSpaceGroup(9, "I 1 1 b", "A 1 1 a | -x+y,-x,z"),

    GeneratedSpaceGroup(10, "P 1 2/m 1", "-x,y,-z; -x,-y,-z", "P 2/m"),
    TransformedSpaceGroup(10, "P 1 1 2/m", "P 1 2/m 1 | y,z,x"),

    GeneratedSpaceGroup(11, "P 1 21/m 1", "-x,y+1/2,-z; -x,-y,-z", "P 21/m"),
    TransformedSpaceGroup(11, "P 1 1 21/m", "P 1 21/m 1 | y,z,x"),

    GeneratedSpaceGroup(12, "C 1 2/m 1", "-x,y,-z; -x,-y,-z", "C 2/m"),
    TransformedSpaceGroup(12, "A 1 2/m 1", "C 1 2/m 1 | -x+z,y,-x", "A 2/m"),
    TransformedSpaceGroup(12, "I 1 2/m 1", "C 1 2/m 1 | -z,y,x-z", "I 2/m"),
    TransformedSpaceGroup(12, "A 1 1 2/m", "C 1 2/m 1 | y,z,x"),
    TransformedSpaceGroup(12, "B 1 1 2/m", "A 1 1 2/m | -y,x-y,z"),
    TransformedSpaceGroup(12, "I 1 1 2/m", "A 1 1 2/m | -x+y,-x,z"),

    GeneratedSpaceGroup(13, "P 1 2/c 1", "-x,y,-z+1/2; -x,-y,-z", "P 2/c"),
    TransformedSpaceGroup(13, "P 1 2/n 1", "P 1 2/c 1 | -x+z,y,-x", "P 2/n"),
    TransformedSpaceGroup(13, "P 1 2/a 1", "P 1 2/c 1 | -z,y,x-z", "P 2/a"),
    TransformedSpaceGroup(13, "P 1 1 2/a", "P 1 2/c 1 | y,z,x"),
    TransformedSpaceGroup(13, "P 1 1 2/n", "P 1 1 2/a | -y,x-y,z"),
    TransformedSpaceGroup(13, "P 1 1 2/b", "P 1 1 2/a | -x+y,-x,z"),

    GeneratedSpaceGroup(14, "P 1 21/c 1", "-x,y+1/2,-z+1/2; -x,-y,-z",
                        "P 21/c"),
    TransformedSpaceGroup(14, "P 1 21/n 1", "P 1 21/c 1 | -x+z,y,-x", "P 21/n"),
    TransformedSpaceGroup(14, "P 1 21/a 1", "P 1 21/c 1 | -z,y,x-z", "P 21/a"),
    TransformedSpaceGroup(14, "P 1 1 21/a", "P 1 21/c 1 | y,z,x"),
    TransformedSpaceGroup(14, "P 1 1 21/n", "P 1 1 21/a | -y,x-y,z"),
    TransformedSpaceGroup(14, "P 1 1 21/b", "P 1 1 21/a | -x+y,-x,z"),

    GeneratedSpaceGroup(15, "C 1 2/c 1", "-x,y,-z+1/2; -x,-y,-z", "C 2/c"),
    TransformedSpaceGroup(15, "A 1 2/n 1", "C 1 2/c 1 | -x+z,y,-x", "A 2/n"),
    TransformedSpaceGroup(15, "I 1 2/a 1", "C 1 2/c 1 | -z,y,x-z", "I 2/a"),
    TransformedSpaceGroup(15, "A 1 1 2/a", "C 1 2/c 1 | y,z,x"),
    TransformedSpaceGroup(15, "B 1 1 2/n", "A 1 1 2/a | -y,x-y,z"),
    TransformedSpaceGroup(15, "I 1 1 2/b", "A 1 1 2/a | -x+y,-x,z")));

/* Orthorhombic space groups have a special subscribe-method in the
 * factory, because each group has potentially 6 different settings.
 * In addition some groups have more than one origin choice.
 */
Subscribe orthorhombic((
    OrthorhombicSpaceGroup(16, "P 2 2 2", "-x,y,-z; x,-y,-z"),
    OrthorhombicSpaceGroup(17, "P 2 2 21", "-x,-y,z+1/2; -x,y,-z+1/2"),
    OrthorhombicSpaceGroup(18, "P 21 21 2", "-x,-y,z; -x+1/2,y+1/2,-z"),
    OrthorhombicSpaceGroup(19, "P 21 21 21",
                           "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2"),
    OrthorhombicSpaceGroup(20, "C 2 2 21", "-x,-y,z+1/2; -x,y,-z+1/2"),
    OrthorhombicSpaceGroup(21, "C 2 2 2", "-x,-y,z; -x,y,-z"),
    OrthorhombicSpaceGroup(22, "F 2 2 2", "-x,-y,z; -x,y,-z"),
    OrthorhombicSpaceGroup(23, "I 2 2 2", "-x,-y,z; -x,y,-z"),
    OrthorhombicSpaceGroup(24, "I 21 21 21",
                           "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2"),
    OrthorhombicSpaceGroup(25, "P m m 2", "-x,-y,z; x,-y,z"),
    OrthorhombicSpaceGroup(26, "P m c 21", "-x,-y,z+1/2; x,-y,z+1/2"),
    OrthorhombicSpaceGroup(27, "P c c 2", "-x,-y,z; x,-y,z+1/2"),
    OrthorhombicSpaceGroup(28, "P m a 2", "-x,-y,z; x+1/2,-y,z"),
    OrthorhombicSpaceGroup(29, "P c a 21", "-x,-y,z+1/2; x+1/2,-y,z"),
    OrthorhombicSpaceGroup(30, "P n c 2", "-x,-y,z; x,-y+1/2,z+1/2"),
    OrthorhombicSpaceGroup(31, "P m n 21", "-x+1/2,-y,z+1/2; x+1/2,-y,z+1/2"),
    OrthorhombicSpaceGroup(32, "P b a 2", "-x,-y,z; x+1/2,-y+1/2,z"),
    OrthorhombicSpaceGroup(33, "P n a 21", "-x,-y,z+1/2; x+1/2,-y+1/2,z"),
    OrthorhombicSpaceGroup(34, "P n n 2", "-x,-y,z; x+1/2,-y+1/2,z+1/2"),
    OrthorhombicSpaceGroup(35, "C m m 2", "-x,-y,z; x,-y,z"),
    OrthorhombicSpaceGroup(36, "C m c 21", "-x,-y,z+1/2; x,-y,z+1/2"),
    OrthorhombicSpaceGroup(37, "C c c 2", "-x,-y,z; x,-y,z+1/2"),
    OrthorhombicSpaceGroup(38, "A m m 2", "-x,-y,z; x,-y,z"),
    OrthorhombicSpaceGroup(39, "A e m 2", "-x,-y,z; x,-y,z+1/2"),
    OrthorhombicSpaceGroup(40, "A m a 2", "-x,-y,z; x+1/2,-y,z"),
    OrthorhombicSpaceGroup(41, "A e a 2", "-x,-y,z; x+1/2,-y+1/2,z"),
    OrthorhombicSpaceGroup(42, "F m m 2", "-x,-y,z; x,-y,z"),
    OrthorhombicSpaceGroup(43, "F d d 2", "-x,-y,z; x+1/4,-y+1/4,z+1/4"),
    OrthorhombicSpaceGroup(44, "I m m 2", "-x,-y,z; x,-y,z"),
    OrthorhombicSpaceGroup(45, "I b a 2", "-x,-y,z; x+1/2,-y+1/2,z"),
    OrthorhombicSpaceGroup(46, "I m a 2", "-x,-y,z; x+1/2,-y,z"),
    OrthorhombicSpaceGroup(47, "P m m m", "-x,-y,z; -x,y,-z; -x,-y,-z"),

    OrthorhombicSpaceGroup(48, "P n n n",
                           "-x,-y,z; -x,y,-z; -x+1/2,-y+1/2,-z+1/2"),
    TransformedOrthorhombicSpaceGroup(48, "P n n n :2",
                                      "P n n n | x-1/4,y-1/4,z-1/4"),

    OrthorhombicSpaceGroup(49, "P c c m", "-x,-y,z; -x,y,-z+1/2; -x,-y,-z"),

    OrthorhombicSpaceGroup(50, "P b a n", "-x,-y,z; -x,y,-z; -x+1/2,-y+1/2,-z"),
    TransformedOrthorhombicSpaceGroup(50, "P b a n :2",
                                      "P b a n | x-1/4,y-1/4,z"),

    OrthorhombicSpaceGroup(51, "P m m a", "-x+1/2,-y,z; -x,y,-z; -x,-y,-z"),
    OrthorhombicSpaceGroup(52, "P n n a",
                           "-x+1/2,-y,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(53, "P m n a",
                           "-x+1/2,-y,z+1/2; -x+1/2,y,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(54, "P c c a", "-x+1/2,-y,z; -x,y,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(55, "P b a m", "-x,-y,z; -x+1/2,y+1/2,-z; -x,-y,-z"),
    OrthorhombicSpaceGroup(56, "P c c n",
                           "-x+1/2,-y+1/2,z; -x,y+1/2,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(57, "P b c m",
                           "-x,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(58, "P n n m",
                           "-x,-y,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z"),

    OrthorhombicSpaceGroup(59, "P m m n",
                           "-x,-y,z; -x+1/2,y+1/2,-z; -x+1/2,-y+1/2,-z"),
    TransformedOrthorhombicSpaceGroup(59, "P m m n :2",
                                      "P m m n | x-1/4,y-1/4,z"),

    OrthorhombicSpaceGroup(60, "P b c n",
                           "-x+1/2,-y+1/2,z+1/2; -x,y,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(61, "P b c a",
                           "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(62, "P n m a",
                           "-x+1/2,-y,z+1/2; -x,y+1/2,-z; -x,-y,-z"),
    OrthorhombicSpaceGroup(63, "C m c m", "-x,-y,z+1/2; -x,y,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(64, "C m c e",
                           "-x,-y+1/2,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(65, "C m m m", "-x,-y,z; -x,y,-z; -x,-y,-z"),
    OrthorhombicSpaceGroup(66, "C c c m", "-x,-y,z; -x,y,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(67, "C m m e", "-x,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z"),
    OrthorhombicSpaceGroup(68, "C c c e",
                           "-x+1/2,-y+1/2,z; -x,y,-z; -x,-y+1/2,-z+1/2"),
    OrthorhombicSpaceGroup(69, "F m m m", "-x,-y,z; -x,y,-z; -x,-y,-z"),

    OrthorhombicSpaceGroup(70, "F d d d",
                           "-x,-y,z; -x,y,-z; -x+1/4,-y+1/4,-z+1/4"),
    TransformedOrthorhombicSpaceGroup(70, "F d d d :2",
                                      "F d d d | x+1/8,y+1/8,z+1/8"),

    OrthorhombicSpaceGroup(71, "I m m m", "-x,-y,z; -x,y,-z; -x,-y,-z"),
    OrthorhombicSpaceGroup(72, "I b a m", "-x,-y,z; -x+1/2,y+1/2,-z; -x,-y,-z"),
    OrthorhombicSpaceGroup(73, "I b c a",
                           "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; -x,-y,-z"),
    OrthorhombicSpaceGroup(74, "I m m a",
                           "-x,-y+1/2,z; -x,y+1/2,-z; -x,-y,-z")));

Subscribe tetragonal((
    GeneratedSpaceGroup(75, "P 4", "-x,-y,z; -y,x,z"),
    GeneratedSpaceGroup(76, "P 41", "-x,-y,z+1/2; -y,x,z+1/4"),
    GeneratedSpaceGroup(77, "P 42", "-x,-y,z; -y,x,z+1/2"),
    GeneratedSpaceGroup(78, "P 43", "-x,-y,z+1/2; -y,x,z+3/4"),
    GeneratedSpaceGroup(79, "I 4", "-x,-y,z; -y,x,z"),
    GeneratedSpaceGroup(80, "I 41", "-x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4"),
    GeneratedSpaceGroup(81, "P -4", "-x,-y,z; y,-x,-z"),
    GeneratedSpaceGroup(82, "I -4", "-x,-y,z; y,-x,-z"),
    GeneratedSpaceGroup(83, "P 4/m", "-x,-y,z; -y,x,z; -x,-y,-z"),
    GeneratedSpaceGroup(84, "P 42/m", "-x,-y,z; -y,x,z+1/2; -x,-y,-z"),

    GeneratedSpaceGroup(85, "P 4/n",
                        "-x,-y,z; -y+1/2,x+1/2,z; -x+1/2,-y+1/2,-z"),
    TransformedSpaceGroup(85, "P 4/n :2", "P 4/n | x+1/4,y-1/4,z"),

    GeneratedSpaceGroup(86, "P 42/n",
                        "-x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,-y+1/2,-z+1/2"),
    TransformedSpaceGroup(86, "P 42/n :2", "P 42/n | x+1/4,y+1/4,z+1/4"),

    GeneratedSpaceGroup(87, "I 4/m", "-x,-y,z; -y,x,z; -x,-y,-z"),

    GeneratedSpaceGroup(
        88, "I 41/a", "-x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; -x,-y+1/2,-z+1/4"),
    TransformedSpaceGroup(88, "I 41/a :2", "I 41/a | x,y+1/4,z+1/8"),

    GeneratedSpaceGroup(89, "P 4 2 2", "-x,-y,z; -y,x,z; -x,y,-z"),
    GeneratedSpaceGroup(90, "P 4 21 2",
                        "-x,-y,z; -y+1/2,x+1/2,z; -x+1/2,y+1/2,-z"),
    GeneratedSpaceGroup(91, "P 41 2 2", "-x,-y,z+1/2; -y,x,z+1/4; -x,y,-z"),
    GeneratedSpaceGroup(92, "P 41 21 2",
                        "-x,-y,z+1/2; -y+1/2,x+1/2,z+1/4; -x+1/2,y+1/2,-z+1/4"),
    GeneratedSpaceGroup(93, "P 42 2 2", "-x,-y,z; -y,x,z+1/2; -x,y,-z"),
    GeneratedSpaceGroup(94, "P 42 21 2",
                        "-x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z+1/2"),
    GeneratedSpaceGroup(95, "P 43 2 2", "-x,-y,z+1/2; -y,x,z+3/4; -x,y,-z"),
    GeneratedSpaceGroup(96, "P 43 21 2",
                        "-x,-y,z+1/2; -y+1/2,x+1/2,z+3/4; -x+1/2,y+1/2,-z+3/4"),
    GeneratedSpaceGroup(97, "I 4 2 2", "-x,-y,z; -y,x,z; -x,y,-z"),
    GeneratedSpaceGroup(98, "I 41 2 2",
                        "-x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; -x+1/2,y,-z+3/4"),
    GeneratedSpaceGroup(99, "P 4 m m", "-x,-y,z; -y,x,z; x,-y,z"),
    GeneratedSpaceGroup(100, "P 4 b m", "-x,-y,z; -y,x,z; x+1/2,-y+1/2,z"),
    GeneratedSpaceGroup(101, "P 42 c m", "-x,-y,z; -y,x,z+1/2; x,-y,z+1/2"),
    GeneratedSpaceGroup(102, "P 42 n m",
                        "-x,-y,z; -y+1/2,x+1/2,z+1/2; x+1/2,-y+1/2,z+1/2"),
    GeneratedSpaceGroup(103, "P 4 c c", "-x,-y,z; -y,x,z; x,-y,z+1/2"),
    GeneratedSpaceGroup(104, "P 4 n c", "-x,-y,z; -y,x,z; x+1/2,-y+1/2,z+1/2"),
    GeneratedSpaceGroup(105, "P 42 m c", "-x,-y,z; -y,x,z+1/2; x,-y,z"),
    GeneratedSpaceGroup(106, "P 42 b c", "-x,-y,z; -y,x,z+1/2; x+1/2,-y+1/2,z"),
    GeneratedSpaceGroup(107, "I 4 m m", "-x,-y,z; -y,x,z; x,-y,z"),
    GeneratedSpaceGroup(108, "I 4 c m", "-x,-y,z; -y,x,z; x,-y,z+1/2"),
    GeneratedSpaceGroup(109, "I 41 m d",
                        "-x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; x,-y,z"),
    GeneratedSpaceGroup(110, "I 41 c d",
                        "-x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; x,-y,z+1/2"),
    GeneratedSpaceGroup(111, "P -4 2 m", "-x,-y,z; y,-x,-z; -x,y,-z"),
    GeneratedSpaceGroup(112, "P -4 2 c", "-x,-y,z; y,-x,-z; -x,y,-z+1/2"),
    GeneratedSpaceGroup(113, "P -4 21 m", "-x,-y,z; y,-x,-z; -x+1/2,y+1/2,-z"),
    GeneratedSpaceGroup(114, "P -4 21 c",
                        "-x,-y,z; y,-x,-z; -x+1/2,y+1/2,-z+1/2"),
    GeneratedSpaceGroup(115, "P -4 m 2", "-x,-y,z; y,-x,-z; x,-y,z"),
    GeneratedSpaceGroup(116, "P -4 c 2", "-x,-y,z; y,-x,-z; x,-y,z+1/2"),
    GeneratedSpaceGroup(117, "P -4 b 2", "-x,-y,z; y,-x,-z; x+1/2,-y+1/2,z"),
    GeneratedSpaceGroup(118, "P -4 n 2",
                        "-x,-y,z; y,-x,-z; x+1/2,-y+1/2,z+1/2"),
    GeneratedSpaceGroup(119, "I -4 m 2", "-x,-y,z; y,-x,-z; x,-y,z"),
    GeneratedSpaceGroup(120, "I -4 c 2", "-x,-y,z; y,-x,-z; x,-y,z+1/2"),
    GeneratedSpaceGroup(121, "I -4 2 m", "-x,-y,z; y,-x,-z; -x,y,-z"),
    GeneratedSpaceGroup(122, "I -4 2 d", "-x,-y,z; y,-x,-z; -x+1/2,y,-z+3/4"),
    GeneratedSpaceGroup(123, "P 4/m m m", "-x,-y,z; -y,x,z; -x,y,-z; -x,-y,-z"),
    GeneratedSpaceGroup(124, "P 4/m c c",
                        "-x,-y,z; -y,x,z; -x,y,-z+1/2; -x,-y,-z"),

    GeneratedSpaceGroup(125, "P 4/n b m",
                        "-x,-y,z; -y,x,z; -x,y,-z; -x+1/2,-y+1/2,-z"),
    TransformedSpaceGroup(125, "P 4/n b m :2", "P 4/n b m | x+1/4,y+1/4,z"),

    GeneratedSpaceGroup(126, "P 4/n n c",
                        "-x,-y,z; -y,x,z; -x,y,-z; -x+1/2,-y+1/2,-z+1/2"),
    TransformedSpaceGroup(126, "P 4/n n c :2", "P 4/n n c | x+1/4,y+1/4,z+1/4"),

    GeneratedSpaceGroup(127, "P 4/m b m",
                        "-x,-y,z; -y,x,z; -x+1/2,y+1/2,-z; -x,-y,-z"),
    GeneratedSpaceGroup(128, "P 4/m n c",
                        "-x,-y,z; -y,x,z; -x+1/2,y+1/2,-z+1/2; -x,-y,-z"),

    GeneratedSpaceGroup(
        129, "P 4/n m m",
        "-x,-y,z; -y+1/2,x+1/2,z; -x+1/2,y+1/2,-z; -x+1/2,-y+1/2,-z"),
    TransformedSpaceGroup(129, "P 4/n m m :2", "P 4/n m m | x+1/4,y-1/4,z"),

    GeneratedSpaceGroup(
        130, "P 4/n c c",
        "-x,-y,z; -y+1/2,x+1/2,z; -x+1/2,y+1/2,-z+1/2; -x+1/2,-y+1/2,-z"),
    TransformedSpaceGroup(130, "P 4/n c c :2", "P 4/n c c | x+1/4,y-1/4,z"),

    GeneratedSpaceGroup(131, "P 42/m m c",
                        "-x,-y,z; -y,x,z+1/2; -x,y,-z; -x,-y,-z"),
    GeneratedSpaceGroup(132, "P 42/m c m",
                        "-x,-y,z; -y,x,z+1/2; -x,y,-z+1/2; -x,-y,-z"),

    GeneratedSpaceGroup(
        133, "P 42/n b c",
        "-x,-y,z; -y+1/2,x+1/2,z+1/2; -x,y,-z+1/2; -x+1/2,-y+1/2,-z+1/2"),
    TransformedSpaceGroup(133, "P 42/n b c :2",
                          "P 42/n b c | x+1/4,y-1/4,z+1/4"),

    GeneratedSpaceGroup(
        134, "P 42/n n m",
        "-x,-y,z; -y+1/2,x+1/2,z+1/2; -x,y,-z; -x+1/2,-y+1/2,-z+1/2"),
    TransformedSpaceGroup(134, "P 42/n n m :2",
                          "P 42/n n m | x+1/4,y-1/4,z+1/4"),

    GeneratedSpaceGroup(135, "P 42/m b c",
                        "-x,-y,z; -y,x,z+1/2; -x+1/2,y+1/2,-z; -x,-y,-z"),
    GeneratedSpaceGroup(
        136, "P 42/m n m",
        "-x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z+1/2; -x,-y,-z"),

    GeneratedSpaceGroup(137, "P 42/n m c",
                        "-x,-y,z; -y+1/2,x+1/2,z+1/2; "
                        "-x+1/2,y+1/2,-z+1/2; -x+1/2,-y+1/2,-z+1/2"),
    TransformedSpaceGroup(137, "P 42/n m c :2",
                          "P 42/n m c | x+1/4,y-1/4,z+1/4"),

    GeneratedSpaceGroup(
        138, "P 42/n c m",
        "-x,-y,z; -y+1/2,x+1/2,z+1/2; -x+1/2,y+1/2,-z; -x+1/2,-y+1/2,-z+1/2"),
    TransformedSpaceGroup(138, "P 42/n c m :2",
                          "P 42/n c m | x+1/4,y-1/4,z+1/4"),

    GeneratedSpaceGroup(139, "I 4/m m m", "-x,-y,z; -y,x,z; -x,y,-z; -x,-y,-z"),
    GeneratedSpaceGroup(140, "I 4/m c m",
                        "-x,-y,z; -y,x,z; -x,y,-z+1/2; -x,-y,-z"),
    GeneratedSpaceGroup(141, "I 41/a m d",
                        "-x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; "
                        "-x+1/2,y,-z+3/4; -x,-y+1/2,-z+1/4"),
    TransformedSpaceGroup(141, "I 41/a m d :2", "I 41/a m d | x,y-1/4,z+1/8"),

    GeneratedSpaceGroup(142, "I 41/a c d",
                        "-x+1/2,-y+1/2,z+1/2; -y,x+1/2,z+1/4; "
                        "-x+1/2,y,-z+1/4; -x,-y+1/2,-z+1/4"),
    TransformedSpaceGroup(142, "I 41/a c d :2", "I 41/a c d | x,y-1/4,z+1/8")));

Subscribe trigonal(
    (GeneratedSpaceGroup(143, "P 3", "-y,x-y,z"),
     GeneratedSpaceGroup(144, "P 31", "-y,x-y,z+1/3"),
     GeneratedSpaceGroup(145, "P 32", "-y,x-y,z+2/3"),

     GeneratedSpaceGroup(146, "R 3", "-y,x-y,z"),
     TransformedSpaceGroup(
         146, "R 3 :r", "R 3 | 2/3x-1/3y-1/3z, 1/3x+1/3y-2/3z, 1/3x+1/3y+1/3z"),

     GeneratedSpaceGroup(147, "P -3", "-y,x-y,z; -x,-y,-z"),

     GeneratedSpaceGroup(148, "R -3", "-y,x-y,z; -x,-y,-z"),
     TransformedSpaceGroup(
         148, "R -3 :r",
         "R -3 | 2/3x-1/3y-1/3z, 1/3x+1/3y-2/3z, 1/3x+1/3y+1/3z"),

     GeneratedSpaceGroup(149, "P 3 1 2", "-y,x-y,z; -y,-x,-z"),
     GeneratedSpaceGroup(150, "P 3 2 1", "-y,x-y,z; y,x,-z"),
     GeneratedSpaceGroup(151, "P 31 1 2", "-y,x-y,z+1/3; -y,-x,-z+2/3"),
     GeneratedSpaceGroup(152, "P 31 2 1", "-y,x-y,z+1/3; y,x,-z"),
     GeneratedSpaceGroup(153, "P 32 1 2", "-y,x-y,z+2/3; -y,-x,-z+1/3"),
     GeneratedSpaceGroup(154, "P 32 2 1", "-y,x-y,z+2/3; y,x,-z"),

     GeneratedSpaceGroup(155, "R 32", "-y,x-y,z; y,x,-z"),
     TransformedSpaceGroup(
         155, "R 32 :r",
         "R 32 | 2/3x-1/3y-1/3z, 1/3x+1/3y-2/3z, 1/3x+1/3y+1/3z"),

     GeneratedSpaceGroup(156, "P 3 m 1", "-y,x-y,z; -y,-x,z"),
     GeneratedSpaceGroup(157, "P 3 1 m", "-y,x-y,z; y,x,z"),
     GeneratedSpaceGroup(158, "P 3 c 1", "-y,x-y,z; -y,-x,z+1/2"),
     GeneratedSpaceGroup(159, "P 3 1 c", "-y,x-y,z; y,x,z+1/2"),

     GeneratedSpaceGroup(160, "R 3 m", "-y,x-y,z; -y,-x,z"),
     TransformedSpaceGroup(
         160, "R 3 m :r",
         "R 3 m | 2/3x-1/3y-1/3z, 1/3x+1/3y-2/3z, 1/3x+1/3y+1/3z"),

     GeneratedSpaceGroup(161, "R 3 c", "-y,x-y,z; -y,-x,z+1/2"),
     TransformedSpaceGroup(
         161, "R 3 c :r",
         "R 3 c | 2/3x-1/3y-1/3z, 1/3x+1/3y-2/3z, 1/3x+1/3y+1/3z"),

     GeneratedSpaceGroup(162, "P -3 1 m", "-y,x-y,z; -y,-x,-z; -x,-y,-z"),
     GeneratedSpaceGroup(163, "P -3 1 c", "-y,x-y,z; -y,-x,-z+1/2; -x,-y,-z"),
     GeneratedSpaceGroup(164, "P -3 m 1", "-y,x-y,z; y,x,-z; -x,-y,-z"),
     GeneratedSpaceGroup(165, "P -3 c 1", "-y,x-y,z; y,x,-z+1/2; -x,-y,-z"),

     GeneratedSpaceGroup(166, "R -3 m", "-y,x-y,z; y,x,-z; -x,-y,-z"),
     TransformedSpaceGroup(
         166, "R -3 m :r",
         "R -3 m | 2/3x-1/3y-1/3z, 1/3x+1/3y-2/3z, 1/3x+1/3y+1/3z"),

     GeneratedSpaceGroup(167, "R -3 c", "-y,x-y,z; y,x,-z+1/2; -x,-y,-z"),
     TransformedSpaceGroup(
         167, "R -3 c :r",
         "R -3 c | 2/3x-1/3y-1/3z, 1/3x+1/3y-2/3z, 1/3x+1/3y+1/3z")));

Subscribe hexagonal((
    GeneratedSpaceGroup(168, "P 6", "-y,x-y,z; -x,-y,z"),
    GeneratedSpaceGroup(169, "P 61", "-y,x-y,z+1/3; -x,-y,z+1/2"),
    GeneratedSpaceGroup(170, "P 65", "-y,x-y,z+2/3; -x,-y,z+1/2"),
    GeneratedSpaceGroup(171, "P 62", "-y,x-y,z+2/3; -x,-y,z"),
    GeneratedSpaceGroup(172, "P 64", "-y,x-y,z+1/3; -x,-y,z"),
    GeneratedSpaceGroup(173, "P 63", "-y,x-y,z; -x,-y,z+1/2"),
    GeneratedSpaceGroup(174, "P -6", "-y,x-y,z; x,y,-z"),
    GeneratedSpaceGroup(175, "P 6/m", "-y,x-y,z; -x,-y,z; -x,-y,-z"),
    GeneratedSpaceGroup(176, "P 63/m", "-y,x-y,z; -x,-y,z+1/2; -x,-y,-z"),
    GeneratedSpaceGroup(177, "P 6 2 2", "-y,x-y,z; -x,-y,z; y,x,-z"),
    GeneratedSpaceGroup(178, "P 61 2 2",
                        "-y,x-y,z+1/3; -x,-y,z+1/2; y,x,-z+1/3"),
    GeneratedSpaceGroup(179, "P 65 2 2",
                        "-y,x-y,z+2/3; -x,-y,z+1/2; y,x,-z+2/3"),
    GeneratedSpaceGroup(180, "P 62 2 2", "-y,x-y,z+2/3; -x,-y,z; y,x,-z+2/3"),
    GeneratedSpaceGroup(181, "P 64 2 2", "-y,x-y,z+1/3; -x,-y,z; y,x,-z+1/3"),
    GeneratedSpaceGroup(182, "P 63 2 2", "-y,x-y,z; -x,-y,z+1/2; y,x,-z"),
    GeneratedSpaceGroup(183, "P 6 m m", "-y,x-y,z; -x,-y,z; -y,-x,z"),
    GeneratedSpaceGroup(184, "P 6 c c", "-y,x-y,z; -x,-y,z; -y,-x,z+1/2"),
    GeneratedSpaceGroup(185, "P 63 c m", "-y,x-y,z; -x,-y,z+1/2; -y,-x,z+1/2"),
    GeneratedSpaceGroup(186, "P 63 m c", "-y,x-y,z; -x,-y,z+1/2; -y,-x,z"),
    GeneratedSpaceGroup(187, "P -6 m 2", "-y,x-y,z; x,y,-z; -y,-x,z"),
    GeneratedSpaceGroup(188, "P -6 c 2", "-y,x-y,z; x,y,-z+1/2; -y,-x,z+1/2"),
    GeneratedSpaceGroup(189, "P -6 2 m", "-y,x-y,z; x,y,-z; y,x,-z"),
    GeneratedSpaceGroup(190, "P -6 2 c", "-y,x-y,z; x,y,-z+1/2; y,x,-z"),
    GeneratedSpaceGroup(191, "P 6/m m m", "-y,x-y,z; x,y,-z; y,x,-z; -x,-y,-z"),
    GeneratedSpaceGroup(192, "P 6/m c c",
                        "-y,x-y,z; -x,-y,z; y,x,-z+1/2; -x,-y,-z"),
    GeneratedSpaceGroup(193, "P 63/m c m",
                        "-y,x-y,z; -x,-y,z+1/2; y,x,-z+1/2; -x,-y,-z"),
    GeneratedSpaceGroup(194, "P 63/m m c",
                        "-y,x-y,z; -x,-y,z+1/2; y,x,-z; -x,-y,-z")));
Subscribe cubic(
    (GeneratedSpaceGroup(195, "P 2 3", "-x,-y,z; -x,y,-z; z,x,y"),
     GeneratedSpaceGroup(196, "F 2 3", "-x,-y,z; -x,y,-z; z,x,y"),
     GeneratedSpaceGroup(197, "I 2 3", "-x,-y,z; -x,y,-z; z,x,y"),
     GeneratedSpaceGroup(198, "P 21 3",
                         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y"),
     GeneratedSpaceGroup(199, "I 21 3",
                         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y"),
     GeneratedSpaceGroup(200, "P m -3", "-x,-y,z; -x,y,-z; z,x,y; -x,-y,-z"),

     GeneratedSpaceGroup(201, "P n -3",
                         "-x,-y,z; -x,y,-z; z,x,y; -x+1/2,-y+1/2,-z+1/2"),
     TransformedSpaceGroup(201, "P n -3 :2", "P n -3 | x+1/4,y+1/4,z+1/4"),

     GeneratedSpaceGroup(202, "F m -3", "-x,-y,z; -x,y,-z; z,x,y; -x,-y,-z"),

     GeneratedSpaceGroup(203, "F d -3",
                         "-x,-y,z; -x,y,-z; z,x,y; -x+1/4,-y+1/4,-z+1/4"),
     TransformedSpaceGroup(203, "F d -3 :2", "F d -3 | x+1/8,y+1/8,z+1/8"),

     GeneratedSpaceGroup(204, "I m -3", "-x,-y,z; -x,y,-z; z,x,y; -x,-y,-z"),
     GeneratedSpaceGroup(205, "P a -3",
                         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; -x,-y,-z"),
     GeneratedSpaceGroup(206, "I a -3",
                         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; -x,-y,-z"),
     GeneratedSpaceGroup(207, "P 4 3 2", "-x,-y,z; -x,y,-z; z,x,y; y,x,-z"),
     GeneratedSpaceGroup(208, "P 42 3 2",
                         "-x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2"),
     GeneratedSpaceGroup(209, "F 4 3 2", "-x,-y,z; -x,y,-z; z,x,y; y,x,-z"),

     GeneratedSpaceGroup(
         210, "F 41 3 2",
         "-x,-y+1/2,z+1/2; -x+1/2,y+1/2,-z; z,x,y; y+3/4,x+1/4,-z+3/4"),
     GeneratedSpaceGroup(211, "I 4 3 2", "-x,-y,z; -x,y,-z; z,x,y; y,x,-z"),
     GeneratedSpaceGroup(
         212, "P 43 3 2",
         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+1/4,x+3/4,-z+3/4"),
     GeneratedSpaceGroup(
         213, "P 41 3 2",
         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4"),
     GeneratedSpaceGroup(
         214, "I 41 3 2",
         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+3/4,x+1/4,-z+1/4"),
     GeneratedSpaceGroup(215, "P -4 3 m", "-x,-y,z; -x,y,-z; z,x,y; y,x,z"),
     GeneratedSpaceGroup(216, "F -4 3 m", "-x,-y,z; -x,y,-z; z,x,y; y,x,z"),
     GeneratedSpaceGroup(217, "I -4 3 m", "-x,-y,z; -x,y,-z; z,x,y; y,x,z"),
     GeneratedSpaceGroup(218, "P -4 3 n",
                         "-x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,z+1/2"),
     GeneratedSpaceGroup(219, "F -4 3 c",
                         "-x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,z+1/2"),
     GeneratedSpaceGroup(
         220, "I -4 3 d",
         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; y+1/4,x+1/4,z+1/4"),
     GeneratedSpaceGroup(221, "P m -3 m",
                         "-x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z"),

     GeneratedSpaceGroup(
         222, "P n -3 n",
         "-x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x+1/2,-y+1/2,-z+1/2"),
     TransformedSpaceGroup(222, "P n -3 n :2", "P n -3 n | x+1/4,y+1/4,z+1/4"),

     GeneratedSpaceGroup(
         223, "P m -3 n",
         "-x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; -x,-y,-z"),

     GeneratedSpaceGroup(224, "P n -3 m",
                         "-x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; "
                         "-x+1/2,-y+1/2,-z+1/2"),
     TransformedSpaceGroup(224, "P n -3 m :2", "P n -3 m | x+1/4,y+1/4,z+1/4"),

     GeneratedSpaceGroup(225, "F m -3 m",
                         "-x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z"),
     GeneratedSpaceGroup(
         226, "F m -3 c",
         "-x,-y,z; -x,y,-z; z,x,y; y+1/2,x+1/2,-z+1/2; -x,-y,-z"),

     GeneratedSpaceGroup(227, "F d -3 m",
                         "-x,-y+1/2,z+1/2; -x+1/2,y+1/2,-z; z,x,y; "
                         "y+3/4,x+1/4,-z+3/4; -x+1/4,-y+1/4,-z+1/4"),
     TransformedSpaceGroup(227, "F d -3 m :2", "F d -3 m | x+1/8,y+1/8,z+1/8"),

     GeneratedSpaceGroup(228, "F d -3 c",
                         "-x,-y+1/2,z+1/2; -x+1/2,y+1/2,-z; z,x,y; "
                         "y+3/4,x+1/4,-z+3/4; -x+3/4,-y+3/4,-z+3/4"),
     TransformedSpaceGroup(228, "F d -3 c :2", "F d -3 c | x+3/8,y+3/8,z+3/8"),

     GeneratedSpaceGroup(229, "I m -3 m",
                         "-x,-y,z; -x,y,-z; z,x,y; y,x,-z; -x,-y,-z"),
     GeneratedSpaceGroup(230, "I a -3 d",
                         "-x+1/2,-y,z+1/2; -x,y+1/2,-z+1/2; z,x,y; "
                         "y+3/4,x+1/4,-z+1/4; -x,-y,-z")));
}
}
}
