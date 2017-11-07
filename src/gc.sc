import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.gatk.HaplotypeCaller
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode
import org.broadinstitute.gatk.utils.commandline.Output
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType

class GvcfWholeGenome extends QScript {
  @Input(doc="The input bam file", shortName = "I")
  var bamFile: List[File] = Nil

  @Input(doc="The reference", shortName="R")
  var referenceFile: File = _

  @Input(doc="dbSNP vcf", shortName="D")
  var dbSNPFile: File = _

  @Output(doc="output", shortName="O")
  var outFile: File = _

  def script(): Unit = {
    val genotyper = new HaplotypeCaller
    genotyper.scatterCount = 200
    genotyper.input_file = this.bamFile
    genotyper.reference_sequence = this.referenceFile
    genotyper.dbsnp = this.dbSNPFile
    genotyper.out = outFile
    genotyper.stand_call_conf = 30
    genotyper.emitRefConfidence = ReferenceConfidenceMode.GVCF
    genotyper.variant_index_type = GATKVCFIndexType.LINEAR
    genotyper.variant_index_parameter = 128000
    add(genotyper)

  }

}